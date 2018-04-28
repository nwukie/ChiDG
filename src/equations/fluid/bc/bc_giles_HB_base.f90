module bc_giles_HB_base
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, FOUR, HALF, ME, CYLINDRICAL,    &
                                      XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, PI
    use mod_fluid,              only: Rgas, cp, gam
    use mod_interpolation,      only: interpolate_linear, interpolate_linear_ad
    use mod_gridspace,          only: linspace
    use mod_dft,                only: dft, idft_eval
    use mod_chimera,            only: find_gq_donor, find_gq_donor_parallel

    use type_point,             only: point_t
    use type_mesh,              only: mesh_t
    use type_bc_state,          only: bc_state_t
    use type_bc_patch,          only: bc_patch_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_face_info,         only: face_info_t, face_info_constructor
    use type_element_info,      only: element_info_t
    use mod_chidg_mpi,          only: IRANK
    use mod_interpolate,        only: interpolate_face_autodiff
    use mpi_f08,                only: MPI_REAL8, MPI_AllReduce, mpi_comm, MPI_INTEGER, MPI_BCast, MPI_MIN, MPI_MAX
    use ieee_arithmetic,        only: ieee_is_nan
    use DNAD_D
    implicit none



    !>  Name: Outlet - 3D Giles
    !!
    !!  Options:
    !!      : Average Pressure
    !!
    !!  Behavior:
    !!      
    !!  References:
    !!              
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2018
    !!
    !---------------------------------------------------------------------------------
    type, public, abstract, extends(bc_state_t) :: giles_HB_base_t

        integer(ik) :: nr = 10
        integer(ik) :: nfourier_space = 14

        real(rk),   allocatable :: r(:)
        real(rk),   allocatable :: theta(:,:)   ! (nr,ntheta)
        real(rk)                :: theta_ref    

        type(element_info_t),   allocatable :: donor(:,:)
        real(rk),               allocatable :: donor_node(:,:,:)

    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: init_bc_coupling     ! Initialize global coupling
        procedure   :: init_bc_postcomm     ! Implement specialized initialization

        procedure   :: get_q_interior

        procedure   :: compute_temporal_dft
        procedure   :: compute_spatial_dft
        procedure   :: compute_absorbing_inlet
        procedure   :: compute_absorbing_outlet
        procedure   :: primitive_to_characteristics
        procedure   :: characteristics_to_primitive
        procedure   :: compute_boundary_average
        procedure   :: analyze_bc_geometry
        procedure   :: initialize_fourier_discretization

    end type giles_HB_base_t
    !*********************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. average_pressure 
    !!  @date   2/8/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(giles_HB_base_t),   intent(inout) :: self

        ! Set name, family
        call self%set_name('Outlet - Giles Quasi3D Unsteady HB')
        call self%set_family('Outlet')

        ! Add functions
        call self%bcproperties%add('Average Pressure',    'Required')
        call self%bcproperties%add('Pitch',               'Required')
        call self%bcproperties%add('Spatial Periodicity', 'Required')

    end subroutine init
    !********************************************************************************




    !>  Initialize boundary group coupling.
    !!
    !!  Call global coupling routine to initialize implicit coupling between each
    !!  element with every other element on the boundary, a result of averaging
    !!  operations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/18/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine init_bc_coupling(self,mesh,group_ID,bc_comm)
        class(giles_HB_base_t), intent(inout)   :: self
        type(mesh_t),           intent(inout)   :: mesh
        integer(ik),            intent(in)      :: group_ID
        type(mpi_comm),         intent(in)      :: bc_comm

        call self%init_bc_coupling_global(mesh,group_ID,bc_comm)

    end subroutine init_bc_coupling
    !********************************************************************************




    !>  Default specialized initialization procedure. This is called from the base 
    !!  bc%init procedure and can be overwritten by derived types to implement 
    !!  specialized initiailization details.
    !!
    !!  By default, this routine does nothing. However, a particular bc_state_t could 
    !!  reimplement this routine to perform some specialized initialization calculations 
    !!  during initialization.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2018
    !!
    !-------------------------------------------------------------------------------------
    subroutine init_bc_postcomm(self,mesh,group_ID,bc_comm)
        class(giles_HB_base_t),   intent(inout)   :: self
        type(mesh_t),                                intent(inout)   :: mesh
        integer(ik),                                 intent(in)      :: group_ID
        type(mpi_comm),                              intent(in)      :: bc_comm

        call self%analyze_bc_geometry(mesh,group_ID,bc_comm)

        call self%initialize_fourier_discretization(mesh,group_ID,bc_comm)

    end subroutine init_bc_postcomm
    !*************************************************************************************



    
    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/25/2018
    !!
    !------------------------------------------------------------------------------------
    subroutine get_q_interior(self,worker,bc_comm,    &
                              density,                &
                              vel1,                   &
                              vel2,                   &
                              vel3,                   &
                              pressure)
        class(giles_HB_base_t),         intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        type(mpi_comm),                 intent(in)      :: bc_comm
        type(AD_D),     allocatable,    intent(inout)   :: density(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: vel1(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: vel2(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: vel3(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: pressure(:,:,:)

        type(AD_D), allocatable, dimension(:,:,:) ::  &
            mom1, mom2, mom3, energy

        integer(ik) :: iradius, itime, nradius, ntheta, ntime, ierr

        ! Define Fourier space discretization to determine
        ! number of theta-samples are being taken
        nradius = size(self%r)
        ntheta  = size(self%theta,2)
        ntime   = worker%time_manager%ntime

        ! Allocate storage for discrete time instances
        allocate(density( nradius,ntheta,ntime),  &
                 mom1(    nradius,ntheta,ntime),  &
                 mom2(    nradius,ntheta,ntime),  &
                 mom3(    nradius,ntheta,ntime),  &
                 vel1(    nradius,ntheta,ntime),  &
                 vel2(    nradius,ntheta,ntime),  &
                 vel3(    nradius,ntheta,ntime),  &
                 energy(  nradius,ntheta,ntime),  &
                 pressure(nradius,ntheta,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError


        ! Perform Fourier decomposition at each radial station.
        do iradius = 1,nradius
            do itime = 1,ntime
                ! Interpolate solution to physical_nodes at current radial station: [ntheta]
                density(iradius,:,itime) = worker%interpolate_field_general('Density',    donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)
                mom1(iradius,:,itime)    = worker%interpolate_field_general('Momentum-1', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)
                mom2(iradius,:,itime)    = worker%interpolate_field_general('Momentum-2', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)
                mom3(iradius,:,itime)    = worker%interpolate_field_general('Momentum-3', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)
                energy(iradius,:,itime)  = worker%interpolate_field_general('Energy',     donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)

                if (worker%coordinate_system() == 'Cylindrical') then
                    mom2(iradius,:,itime) = mom2(iradius,:,itime)/self%r(iradius)  ! convert to tangential momentum
                end if
            end do
        end do

        ! Compute velocities and pressure at each time
        vel1 = mom1/density
        vel2 = mom2/density
        vel3 = mom3/density
        pressure = (gam-ONE)*(energy - HALF*( mom1*mom1 + mom2*mom2 + mom3*mom3 )/density )

    end subroutine get_q_interior
    !************************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/12/2018
    !!
    !------------------------------------------------------------------------------------
    subroutine compute_temporal_dft(self,worker,bc_comm,                    &
                                    density_grid, vel1_grid, vel2_grid, vel3_grid, pressure_grid, &
                                    density_Ft_real,  density_Ft_imag,      &
                                    vel1_Ft_real,     vel1_Ft_imag,         &
                                    vel2_Ft_real,     vel2_Ft_imag,         &
                                    vel3_Ft_real,     vel3_Ft_imag,         &
                                    pressure_Ft_real, pressure_Ft_imag)
        class(giles_HB_base_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        type(mpi_comm),                             intent(in)      :: bc_comm
        type(AD_D),     allocatable,                intent(inout)   :: density_grid(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel1_grid(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel2_grid(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel3_grid(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: pressure_grid(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: density_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: density_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel1_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel1_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel2_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel2_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel3_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel3_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: pressure_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: pressure_Ft_imag(:,:,:)


        type(AD_D), allocatable,    dimension(:)    ::                                          &
            density_real_tmp, vel1_real_tmp, vel2_real_tmp, vel3_real_tmp, pressure_real_tmp,   &
            density_imag_tmp, vel1_imag_tmp, vel2_imag_tmp, vel3_imag_tmp, pressure_imag_tmp

        type(AD_D)  :: density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, c_bar

        integer(ik) :: nradius, ntheta, iradius, itheta, imode, itime, ntime, ierr


        ! Define Fourier space discretization to determine
        ! number of theta-samples are being taken
        ntheta  = size(self%theta,2)
        nradius = size(self%r)
        ntime   = worker%mesh%domain(worker%element_info%idomain_l)%elems(worker%element_info%ielement_l)%ntime

        ! Allocate storage for temporal dft
        allocate(density_Ft_real( nradius,ntheta,ntime), density_Ft_imag( nradius,ntheta,ntime),  &
                 vel1_Ft_real(    nradius,ntheta,ntime), vel1_Ft_imag(    nradius,ntheta,ntime),  &
                 vel2_Ft_real(    nradius,ntheta,ntime), vel2_Ft_imag(    nradius,ntheta,ntime),  &
                 vel3_Ft_real(    nradius,ntheta,ntime), vel3_Ft_imag(    nradius,ntheta,ntime),  &
                 pressure_Ft_real(nradius,ntheta,ntime), pressure_Ft_imag(nradius,ntheta,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Perform Fourier decomposition at each radial station.
        do iradius = 1,nradius
            do itheta = 1,ntheta
                call dft(density_grid(iradius,itheta,:),  ZERO*density_grid(iradius,itheta,:),  density_real_tmp,  density_imag_tmp )
                call dft(vel1_grid(iradius,itheta,:),     ZERO*vel1_grid(iradius,itheta,:),     vel1_real_tmp,     vel1_imag_tmp    )
                call dft(vel2_grid(iradius,itheta,:),     ZERO*vel2_grid(iradius,itheta,:),     vel2_real_tmp,     vel2_imag_tmp    )
                call dft(vel3_grid(iradius,itheta,:),     ZERO*vel3_grid(iradius,itheta,:),     vel3_real_tmp,     vel3_imag_tmp    )
                call dft(pressure_grid(iradius,itheta,:), ZERO*pressure_grid(iradius,itheta,:), pressure_real_tmp, pressure_imag_tmp)

                density_Ft_real( iradius,itheta,:) = density_real_tmp
                density_Ft_imag( iradius,itheta,:) = density_imag_tmp
                vel1_Ft_real(    iradius,itheta,:) = vel1_real_tmp
                vel1_Ft_imag(    iradius,itheta,:) = vel1_imag_tmp
                vel2_Ft_real(    iradius,itheta,:) = vel2_real_tmp
                vel2_Ft_imag(    iradius,itheta,:) = vel2_imag_tmp
                vel3_Ft_real(    iradius,itheta,:) = vel3_real_tmp
                vel3_Ft_imag(    iradius,itheta,:) = vel3_imag_tmp
                pressure_Ft_real(iradius,itheta,:) = pressure_real_tmp
                pressure_Ft_imag(iradius,itheta,:) = pressure_imag_tmp

            end do !itheta
        end do !iradius

    end subroutine compute_temporal_dft
    !********************************************************************************


!    !>
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   4/12/2018
!    !!
!    !------------------------------------------------------------------------------------
!    subroutine compute_temporal_dft(self,worker,bc_comm,                    &
!                                    density_Ft_real,  density_Ft_imag,      &
!                                    vel1_Ft_real,     vel1_Ft_imag,         &
!                                    vel2_Ft_real,     vel2_Ft_imag,         &
!                                    vel3_Ft_real,     vel3_Ft_imag,         &
!                                    pressure_Ft_real, pressure_Ft_imag)
!        class(giles_HB_base_t),  intent(inout)   :: self
!        type(chidg_worker_t),                       intent(inout)   :: worker
!        type(mpi_comm),                             intent(in)      :: bc_comm
!        type(AD_D),     allocatable,                intent(inout)   :: density_Ft_real(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: density_Ft_imag(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: vel1_Ft_real(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: vel1_Ft_imag(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: vel2_Ft_real(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: vel2_Ft_imag(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: vel3_Ft_real(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: vel3_Ft_imag(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: pressure_Ft_real(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: pressure_Ft_imag(:,:,:)
!
!
!        type(AD_D), allocatable,    dimension(:)    ::                                          &
!            density_real_tmp, vel1_real_tmp, vel2_real_tmp, vel3_real_tmp, pressure_real_tmp,   &
!            density_imag_tmp, vel1_imag_tmp, vel2_imag_tmp, vel3_imag_tmp, pressure_imag_tmp
!
!
!        type(AD_D), allocatable, dimension(:,:)     ::  &
!            density, mom1, mom2, mom3, energy, vel1, vel2, vel3, pressure
!
!        type(AD_D)  :: density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, c_bar
!
!        integer(ik) :: nradius, ntheta, iradius, itheta, imode, itime, ntime, ierr
!
!
!
!        ! Define Fourier space discretization to determine
!        ! number of theta-samples are being taken
!        ntheta  = size(self%theta,2)
!        nradius = size(self%r)
!        ntime   = worker%mesh%domain(worker%element_info%idomain_l)%elems(worker%element_info%ielement_l)%ntime
!
!        ! Allocate storage for discrete time instances
!        allocate(density( ntheta,ntime),  &
!                 mom1(    ntheta,ntime),  &
!                 mom2(    ntheta,ntime),  &
!                 mom3(    ntheta,ntime),  &
!                 vel1(    ntheta,ntime),  &
!                 vel2(    ntheta,ntime),  &
!                 vel3(    ntheta,ntime),  &
!                 energy(  ntheta,ntime),  &
!                 pressure(ntheta,ntime), stat=ierr)
!        if (ierr /= 0) call AllocationError
!        
!        ! Allocate storage for temporal dft
!        allocate(density_Ft_real( nradius,ntheta,ntime), density_Ft_imag( nradius,ntheta,ntime),  &
!                 vel1_Ft_real(    nradius,ntheta,ntime), vel1_Ft_imag(    nradius,ntheta,ntime),  &
!                 vel2_Ft_real(    nradius,ntheta,ntime), vel2_Ft_imag(    nradius,ntheta,ntime),  &
!                 vel3_Ft_real(    nradius,ntheta,ntime), vel3_Ft_imag(    nradius,ntheta,ntime),  &
!                 pressure_Ft_real(nradius,ntheta,ntime), pressure_Ft_imag(nradius,ntheta,ntime), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!        ! Perform Fourier decomposition at each radial station.
!        do iradius = 1,nradius
!            do itime = 1,ntime
!                ! Interpolate solution to physical_nodes at current radial station: [ntheta]
!                density(:,itime) = worker%interpolate_field_general('Density',    donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)
!                mom1(:,itime)    = worker%interpolate_field_general('Momentum-1', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)
!                mom2(:,itime)    = worker%interpolate_field_general('Momentum-2', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)
!                mom3(:,itime)    = worker%interpolate_field_general('Momentum-3', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)
!                energy(:,itime)  = worker%interpolate_field_general('Energy',     donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)
!
!                if (worker%coordinate_system() == 'Cylindrical') then
!                    mom2(:,itime) = mom2(:,itime)/self%r(iradius)  ! convert to tangential momentum
!                end if
!            end do
!
!            ! Compute velocities and pressure at each time
!            vel1 = mom1/density
!            vel2 = mom2/density
!            vel3 = mom3/density
!            pressure = (gam-ONE)*(energy - HALF*( mom1*mom1 + mom2*mom2 + mom3*mom3 )/density )
!
!
!            ! Temporal dft at each theta location
!            do itheta = 1,ntheta
!                call dft(density(itheta,:),  ZERO*density(itheta,:),  density_real_tmp,  density_imag_tmp )
!                call dft(vel1(itheta,:),     ZERO*vel1(itheta,:),     vel1_real_tmp,     vel1_imag_tmp    )
!                call dft(vel2(itheta,:),     ZERO*vel2(itheta,:),     vel2_real_tmp,     vel2_imag_tmp    )
!                call dft(vel3(itheta,:),     ZERO*vel3(itheta,:),     vel3_real_tmp,     vel3_imag_tmp    )
!                call dft(pressure(itheta,:), ZERO*pressure(itheta,:), pressure_real_tmp, pressure_imag_tmp)
!
!                density_Ft_real( iradius,itheta,:) = density_real_tmp
!                density_Ft_imag( iradius,itheta,:) = density_imag_tmp
!                vel1_Ft_real(    iradius,itheta,:) = vel1_real_tmp
!                vel1_Ft_imag(    iradius,itheta,:) = vel1_imag_tmp
!                vel2_Ft_real(    iradius,itheta,:) = vel2_real_tmp
!                vel2_Ft_imag(    iradius,itheta,:) = vel2_imag_tmp
!                vel3_Ft_real(    iradius,itheta,:) = vel3_real_tmp
!                vel3_Ft_imag(    iradius,itheta,:) = vel3_imag_tmp
!                pressure_Ft_real(iradius,itheta,:) = pressure_real_tmp
!                pressure_Ft_imag(iradius,itheta,:) = pressure_imag_tmp
!
!            end do !itheta
!
!        end do !iradius
!
!    end subroutine compute_temporal_dft
!    !********************************************************************************

!    !>
!    !!
!    !!
!    !!
!    !!
!    !!
!    !---------------------------------------------------------------------------
!    subroutine compute_temporal_idft(self,worker,bc_comm,   &
!                                     density
!
!        ! Inverse DFT of temporal Fourier modes to give primitive variables
!        ! at quarature nodes for the current time instance.
!        density_bc_tmp  = [ZERO*density_Fts_real_gq(1,1,1)]
!        vel1_bc_tmp     = [ZERO*density_Fts_real_gq(1,1,1)]
!        vel2_bc_tmp     = [ZERO*density_Fts_real_gq(1,1,1)]
!        vel3_bc_tmp     = [ZERO*density_Fts_real_gq(1,1,1)]
!        pressure_bc_tmp = [ZERO*density_Fts_real_gq(1,1,1)]
!        do igq = 1,size(coords)
!            ! **** WARNING: probably want ipdft_eval here ****
!            call idft_eval(density_t_real(igq,:),   &
!                           density_t_imag(igq,:),   &
!                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
!                           density_bc_tmp,          &
!                           expect_zero)
!
!            call idft_eval(vel1_t_real(igq,:),      &
!                           vel1_t_imag(igq,:),      &
!                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
!                           vel1_bc_tmp,        &
!                           expect_zero)
!
!            call idft_eval(vel2_t_real(igq,:),      &
!                           vel2_t_imag(igq,:),      &
!                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
!                           vel2_bc_tmp,        &
!                           expect_zero)
!
!            call idft_eval(vel3_t_real(igq,:),      &
!                           vel3_t_imag(igq,:),      &
!                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
!                           vel3_bc_tmp,        &
!                           expect_zero)
!
!            call idft_eval(pressure_t_real(igq,:),  &
!                           pressure_t_imag(igq,:),  &
!                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
!                           pressure_bc_tmp,    &
!                           expect_zero)
!
!            ! Accumulate contribution from unsteady modes
!            density_bc(igq)  = density_bc(igq)  + density_bc_tmp(1)
!            vel1_bc(igq)     = vel1_bc(igq)     + vel1_bc_tmp(1)
!            vel2_bc(igq)     = vel2_bc(igq)     + vel2_bc_tmp(1)
!            vel3_bc(igq)     = vel3_bc(igq)     + vel3_bc_tmp(1)
!            pressure_bc(igq) = pressure_bc(igq) + pressure_bc_tmp(1)
!        end do
!        
!    end subroutine compute_temporal_idft



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_spatial_dft(self,worker,bc_comm,                   &
                                   density_Ft_real,  density_Ft_imag,     &
                                   vel1_Ft_real,     vel1_Ft_imag,        &
                                   vel2_Ft_real,     vel2_Ft_imag,        &
                                   vel3_Ft_real,     vel3_Ft_imag,        &
                                   pressure_Ft_real, pressure_Ft_imag,    &
                                   density_Fts_real,  density_Fts_imag,   &
                                   vel1_Fts_real,     vel1_Fts_imag,      &
                                   vel2_Fts_real,     vel2_Fts_imag,      &
                                   vel3_Fts_real,     vel3_Fts_imag,      &
                                   pressure_Fts_real, pressure_Fts_imag)
        class(giles_HB_base_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        type(mpi_comm),                             intent(in)      :: bc_comm
        type(AD_D),     allocatable,                intent(inout)   :: density_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: density_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel1_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel1_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel2_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel2_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel3_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel3_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: pressure_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: pressure_Ft_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: density_Fts_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: density_Fts_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel1_Fts_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel1_Fts_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel2_Fts_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel2_Fts_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel3_Fts_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel3_Fts_imag(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: pressure_Fts_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: pressure_Fts_imag(:,:,:)

        type(AD_D), allocatable,    dimension(:)    ::                                          &
            density_real_tmp, vel1_real_tmp, vel2_real_tmp, vel3_real_tmp, pressure_real_tmp,   &
            density_imag_tmp, vel1_imag_tmp, vel2_imag_tmp, vel3_imag_tmp, pressure_imag_tmp

        integer(ik)             :: nradius, ntheta, iradius, itheta, imode, itime, ntime, ierr
        real(rk)                :: shift_r, shift_i
        real(rk),   allocatable :: spatial_periodicity(:)

        ! Define Fourier space discretization to determine
        ! number of theta-samples being taken
        ntheta  = size(self%theta,2)
        nradius = size(self%r)
        ntime   = worker%mesh%domain(worker%element_info%idomain_l)%elems(worker%element_info%ielement_l)%ntime
        spatial_periodicity = self%bcproperties%compute('Spatial Periodicity', time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])
        
        ! Allocate storage in result
        allocate(density_Fts_real( nradius,ntheta,ntime), density_Fts_imag( nradius,ntheta,ntime),  &
                 vel1_Fts_real(    nradius,ntheta,ntime), vel1_Fts_imag(    nradius,ntheta,ntime),  &
                 vel2_Fts_real(    nradius,ntheta,ntime), vel2_Fts_imag(    nradius,ntheta,ntime),  &
                 vel3_Fts_real(    nradius,ntheta,ntime), vel3_Fts_imag(    nradius,ntheta,ntime),  &
                 pressure_Fts_real(nradius,ntheta,ntime), pressure_Fts_imag(nradius,ntheta,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Perform Fourier decomposition at each radial station.
        do iradius = 1,nradius
            do itime = 1,ntime

                ! DFT in space
                call dft(density_Ft_real( iradius,:,itime), density_Ft_imag( iradius,:,itime), density_real_tmp,  density_imag_tmp )
                call dft(vel1_Ft_real(    iradius,:,itime), vel1_Ft_imag(    iradius,:,itime), vel1_real_tmp,     vel1_imag_tmp    )
                call dft(vel2_Ft_real(    iradius,:,itime), vel2_Ft_imag(    iradius,:,itime), vel2_real_tmp,     vel2_imag_tmp    )
                call dft(vel3_Ft_real(    iradius,:,itime), vel3_Ft_imag(    iradius,:,itime), vel3_real_tmp,     vel3_imag_tmp    )
                call dft(pressure_Ft_real(iradius,:,itime), pressure_Ft_imag(iradius,:,itime), pressure_real_tmp, pressure_imag_tmp)

                ! Adjust Fourier coefficients so their phase is relative to self%theta_ref
                ! instead of the minimum theta of the transform.
                !
                !       q(relative to theta_ref) = q(relative to theta_min) * e^(j 2pi imode delta_theta/spatial_periodicity)
                !
                ! NOTE: self%theta(:,1) are defined to be the DFT-theta_min at each radius
                !
                do imode = 1,size(density_real_tmp)
                    shift_r = realpart(exp(cmplx(ZERO,ONE)*real(imode-1,rk)*TWO*PI*(self%theta_ref-self%theta(iradius,1))/spatial_periodicity(1)))
                    shift_i = imagpart(exp(cmplx(ZERO,ONE)*real(imode-1,rk)*TWO*PI*(self%theta_ref-self%theta(iradius,1))/spatial_periodicity(1)))

                    density_Fts_real( iradius,imode,itime) = density_real_tmp(imode)*shift_r  - density_imag_tmp(imode)*shift_i
                    vel1_Fts_real(    iradius,imode,itime) = vel1_real_tmp(imode)*shift_r     - vel1_imag_tmp(imode)*shift_i
                    vel2_Fts_real(    iradius,imode,itime) = vel2_real_tmp(imode)*shift_r     - vel2_imag_tmp(imode)*shift_i
                    vel3_Fts_real(    iradius,imode,itime) = vel3_real_tmp(imode)*shift_r     - vel3_imag_tmp(imode)*shift_i
                    pressure_Fts_real(iradius,imode,itime) = pressure_real_tmp(imode)*shift_r - pressure_imag_tmp(imode)*shift_i

                    density_Fts_imag( iradius,imode,itime) = density_imag_tmp(imode)*shift_r  + density_real_tmp(imode)*shift_i
                    vel1_Fts_imag(    iradius,imode,itime) = vel1_imag_tmp(imode)*shift_r     + vel1_real_tmp(imode)*shift_i
                    vel2_Fts_imag(    iradius,imode,itime) = vel2_imag_tmp(imode)*shift_r     + vel2_real_tmp(imode)*shift_i
                    vel3_Fts_imag(    iradius,imode,itime) = vel3_imag_tmp(imode)*shift_r     + vel3_real_tmp(imode)*shift_i
                    pressure_Fts_imag(iradius,imode,itime) = pressure_imag_tmp(imode)*shift_r + pressure_real_tmp(imode)*shift_i
                end do !imode

            end do !itime
        end do !iradius

        print*, 'WARNING: need scaling since dft is only over a single passage.'
        print*, 'WARNING: check correct pitch in phase shift.'

    end subroutine compute_spatial_dft
    !*********************************************************************************



!    !>  GILES' FORMULATION
!    !!
!    !!
!    !!
!    !!
!    !--------------------------------------------------------------------------------
!    subroutine compute_absorbing_inlet(self,worker,bc_comm, &
!                                       density_real,    &
!                                       density_imag,    &
!                                       vel1_real,       &
!                                       vel1_imag,       &
!                                       vel2_real,       &
!                                       vel2_imag,       &
!                                       vel3_real,       &
!                                       vel3_imag,       &
!                                       pressure_real,   &
!                                       pressure_imag)
!        class(giles_HB_base_t),  intent(inout)   :: self
!        type(chidg_worker_t),    intent(inout)   :: worker
!        type(mpi_comm),          intent(in)      :: bc_comm
!        type(AD_D),              intent(inout)   :: density_real(:,:,:)
!        type(AD_D),              intent(inout)   :: density_imag(:,:,:)
!        type(AD_D),              intent(inout)   :: vel1_real(:,:,:)
!        type(AD_D),              intent(inout)   :: vel1_imag(:,:,:)
!        type(AD_D),              intent(inout)   :: vel2_real(:,:,:)
!        type(AD_D),              intent(inout)   :: vel2_imag(:,:,:)
!        type(AD_D),              intent(inout)   :: vel3_real(:,:,:)
!        type(AD_D),              intent(inout)   :: vel3_imag(:,:,:)
!        type(AD_D),              intent(inout)   :: pressure_real(:,:,:)
!        type(AD_D),              intent(inout)   :: pressure_imag(:,:,:)
!
!
!        ! Storage at quadrature nodes
!        type(AD_D), allocatable, dimension(:)   ::                                      &
!            c1,    c2,    c3,    c4,    c5,                                             &
!            c1_3d, c2_3d, c3_3d, c4_3d, c5_3d,                                          &
!            c1_1d, c2_1d, c3_1d, c4_1d, c5_1d,                                          &
!            c_bar, ddensity, dvel1, dvel2, dvel3, dpressure, expect_zero,               &
!            density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar
!
!        type(AD_D), allocatable, dimension(:,:,:) ::        &
!            c1_real, c2_real, c3_real, c4_real, c5_real,    &
!            c1_imag, c2_imag, c3_imag, c4_imag, c5_imag
!
!        type(AD_D)  :: pressure_avg, vel1_avg, vel2_avg, vel3_avg, density_avg, c_avg, T_avg,       &
!                       density_bar_r, vel1_bar_r, vel2_bar_r, vel3_bar_r, pressure_bar_r, c_bar_r,  &
!                       A3_real, A3_imag, A4_real, A4_imag, A4_denom, A4_real_num, A4_imag_num,      &
!                       beta, s_real, s_imag, s_arg, lambda, vmag
!
!        real(rk),       allocatable, dimension(:)   :: PT, TT, n1, n2, n3, nmag, pitch
!        real(rk)                                    :: theta_offset, omega, lm, a1, a2, a3, a4
!        integer(ik)                                 :: iradius, igq, ierr, itheta, ntheta, itime, ntime
!
!        pitch  = self%bcproperties%compute('Pitch',worker%time(),worker%coords())
!
!        ! Get spatio-temporal average at radial stations
!        density_bar  = density_real(:,1,1)
!        vel1_bar     = vel1_real(:,1,1)
!        vel2_bar     = vel2_real(:,1,1)
!        vel3_bar     = vel3_real(:,1,1)
!        pressure_bar = pressure_real(:,1,1)
!
!
!        ! Compute Fourier decomposition at set of radial stations: 
!        call self%primitive_to_characteristics(worker,bc_COMM,                  &
!                                               density_real,  density_imag,     &
!                                               vel1_real,     vel1_imag,        &
!                                               vel2_real,     vel2_imag,        &
!                                               vel3_real,     vel3_imag,        &
!                                               pressure_real, pressure_imag,    &
!                                               c1_real,       c1_imag,          &
!                                               c2_real,       c2_imag,          &
!                                               c3_real,       c3_imag,          &
!                                               c4_real,       c4_imag,          &
!                                               c5_real,       c5_imag)
!
!        ! Handle temporal average(steady) nonreflecting condition (:,:,1)
!        ntheta = size(c5_real,2)
!        ntime  = size(c5_real,3)
!        do iradius = 1,size(self%r)
!            ! Get average parts
!            density_bar_r  = density_bar(iradius)
!            vel1_bar_r     = vel1_bar(iradius)
!            vel2_bar_r     = vel2_bar(iradius)
!            vel3_bar_r     = vel3_bar(iradius)
!            pressure_bar_r = pressure_bar(iradius)
!            c_bar_r        = sqrt(gam*pressure_bar_r/density_bar_r)
!
!            ! starting with 2 here because the first mode is treated with 1D characteristics
!            do itheta = 2,ntheta
!                ! Account for sign(mode) in the calculation of beta. The second half of the
!                ! modes are negative frequencies.
!                if (itheta <= (ntheta-1)/2 + 1) then
!                    beta = sqrt(c_bar_r*c_bar_r  -  (vel3_bar_r*vel3_bar_r + vel2_bar_r*vel2_bar_r))
!                else if (itheta > (ntheta-1)/2 + 1) then
!                    beta = -sqrt(c_bar_r*c_bar_r  -  (vel3_bar_r*vel3_bar_r + vel2_bar_r*vel2_bar_r))
!                end if
!
!                ! The imaginary part of beta has already been accounted for in
!                ! the expressions for A2 and A3
!                A3_real = -TWO*vel3_bar_r*vel2_bar_r/(vel2_bar_r*vel2_bar_r + beta*beta)
!                A3_imag = -TWO*beta*vel3_bar_r/(vel2_bar_r*vel2_bar_r + beta*beta)
!
!                A4_real = (beta*beta - vel2_bar_r*vel2_bar_r)/(beta*beta + vel2_bar_r*vel2_bar_r)
!                A4_imag = -TWO*beta*vel2_bar_r/(beta*beta + vel2_bar_r*vel2_bar_r)
!
!                c1_real(iradius,itheta,1) = ZERO
!                c1_imag(iradius,itheta,1) = ZERO
!
!                c2_real(iradius,itheta,1) = ZERO
!                c2_imag(iradius,itheta,1) = ZERO
!
!                c3_real(iradius,itheta,1) = -(ONE/(c_bar_r + vel3_bar_r))*(vel2_bar_r*c5_real(iradius,itheta,1) - beta*c5_imag(iradius,itheta,1))
!                c3_imag(iradius,itheta,1) = -(ONE/(c_bar_r + vel3_bar_r))*(vel2_bar_r*c5_imag(iradius,itheta,1) + beta*c5_real(iradius,itheta,1))
!
!                c4_real(iradius,itheta,1) = (ONE/((c_bar_r+vel3_bar_r)**TWO))*((vel2_bar_r*vel2_bar_r - beta*beta)*c5_real(iradius,itheta,1) - TWO*beta*c5_imag(iradius,itheta,1))
!                c4_imag(iradius,itheta,1) = (ONE/((c_bar_r+vel3_bar_r)**TWO))*((vel2_bar_r*vel2_bar_r - beta*beta)*c5_imag(iradius,itheta,1) + TWO*beta*c5_real(iradius,itheta,1))
!
!            end do !itheta
!        end do !iradius
!
!
!        ! Handle unsteady nonreflecting condition (:,:,2:)
!        print*, 'WARNING! HARDCODED ACCESS TO FIRST FREQUENCY ONLY!!!'
!        print*, 'WARNING! CHECK DEFINITION OF lm PITCH!'
!        do iradius = 1,size(self%r)
!            ! Get average parts
!            density_bar_r  = density_bar(iradius)
!            vel1_bar_r     = vel1_bar(iradius)
!            vel2_bar_r     = vel2_bar(iradius)
!            vel3_bar_r     = vel3_bar(iradius)
!            pressure_bar_r = pressure_bar(iradius)
!            c_bar_r        = sqrt(gam*pressure_bar_r/density_bar_r)
!
!
!            ! starting with 2 here because the first mode is treated with 1D characteristics
!            do itheta = 1,ntheta
!                do itime = 2,ntime
!
!                    omega = worker%time_manager%freqs(1)
!                    lm = TWO*PI*real(itheta-1,rk)/pitch(1)
!                    lambda = lm*c_bar_r/omega
!                    s_arg = ONE - (c_bar_r*c_bar_r - vel3_bar_r*vel3_bar_r)*lambda*lambda/((c_bar_r - vel2_bar_r*lambda)**TWO)
!
!                    if (s_arg >= ZERO) then
!                        s_real = sqrt(s_arg)
!                        s_imag = ZERO*s_arg
!                    else
!                        s_real = ZERO*s_arg
!                        s_imag = sqrt(-s_arg) ! minus because we have already factored out sqrt(-1) = i into the formula below
!                    end if
!
!                    A3_real =  (c_bar_r - vel3_bar_r)*lambda*(ONE + s_real)/((c_bar_r-vel2_bar_r*lambda)*((ONE+s_real)**TWO + s_imag*s_imag))
!                    A3_imag = -(c_bar_r - vel3_bar_r)*lambda*s_imag/((c_bar_r-vel2_bar_r*lambda)*((ONE+s_real)**TWO + s_imag*s_imag))
!
!                    !A4_denom = (c_bar_r - vel2_bar_r*lambda)*(((ONE+s_real)**TWO - s_imag*s_imag)**TWO + FOUR*s_imag*s_imag*((ONE+s_real)**TWO))
!                    A4_denom = ((c_bar_r - vel2_bar_r*lambda)**TWO)*(((ONE+s_real)**TWO - s_imag*s_imag)**TWO + FOUR*s_imag*s_imag*((ONE+s_real)**TWO))
!                    A4_real_num =  ((c_bar_r - vel3_bar_r)**TWO)*lambda*lambda*((ONE+s_real)**TWO - s_imag*s_imag)
!                    A4_imag_num = -((c_bar_r - vel3_bar_r)**TWO)*lambda*lambda*TWO*s_imag*(ONE+s_real)
!                    A4_real = A4_real_num/A4_denom
!                    A4_imag = A4_imag_num/A4_denom
!
!                    c1_real(iradius,itheta,itime) = ZERO
!                    c1_imag(iradius,itheta,itime) = ZERO
!
!                    c2_real(iradius,itheta,itime) = ZERO
!                    c2_imag(iradius,itheta,itime) = ZERO
!
!                    c3_real(iradius,itheta,itime) = A3_real*c5_real(iradius,itheta,itime) - A3_imag*c5_imag(iradius,itheta,itime)
!                    c3_imag(iradius,itheta,itime) = A3_real*c5_imag(iradius,itheta,itime) + A3_imag*c5_real(iradius,itheta,itime)
!                    
!                    c4_real(iradius,itheta,itime) = A4_real*c5_real(iradius,itheta,itime) - A4_imag*c5_imag(iradius,itheta,itime)
!                    c4_imag(iradius,itheta,itime) = A4_real*c5_imag(iradius,itheta,itime) + A4_imag*c5_real(iradius,itheta,itime)
!
!                    ! Add incoming amplitude
!                    if (itheta == 2 .and. itime == 2) then
!                        a3 = 0.
!                        a4 = 0.
!
!                        A3_real = c_bar_r*(vel3_bar_r - c_bar_r)*(vel3_bar_r*s_real - c_bar_r)/((c_bar_r-vel2_bar_r*lambda)*((vel3_bar_r*s_real-c_bar_r)**TWO + vel3_bar_r*vel3_bar_r*s_imag*s_imag))
!                        A3_imag = -vel3_bar_r*c_bar_r*s_imag*(vel3_bar_r-c_bar_r)/((c_bar_r-vel2_bar_r*lambda)*((vel3_bar_r*s_real-c_bar_r)**TWO + vel3_bar_r*vel3_bar_r*s_imag*s_imag))
!                        c3_real(iradius,itheta,itime) = c3_real(iradius,itheta,itime) + A3_real*a3
!                        c3_imag(iradius,itheta,itime) = c3_imag(iradius,itheta,itime) + A3_imag*a3
!
!                        A4_denom    = ((c_bar_r - vel2_bar_r*lambda)**TWO)*( (vel3_bar_r*(s_real*s_real-s_imag*s_imag) - s_real*(c_bar_r - vel3_bar_r) - c_bar_r)**TWO + &
!                                                                             (TWO*vel3_bar_r*s_real*s_imag - s_imag*(c_bar_r-vel3_bar_r))**TWO )
!                        A4_real_num = TWO*lambda*vel3_bar_r*c_bar_r*(c_bar_r-vel3_bar_r) * (vel3_bar_r*(s_real*s_real-s_imag*s_imag) - s_real*(c_bar_r - vel3_bar_r) - c_bar_r)
!                        A4_imag_num = TWO*lambda*vel3_bar_r*c_bar_r*(c_bar_r-vel3_bar_r) * (-(TWO*vel3_bar_r*s_real*s_imag - s_imag*(c_bar_r-vel3_bar_r)))
!                        A4_real = A4_real_num/A4_denom
!                        A4_imag = A4_imag_num/A4_denom
!                        c4_real(iradius,itheta,itime) = c4_real(iradius,itheta,itime) + A4_real*a3
!                        c4_imag(iradius,itheta,itime) = c4_imag(iradius,itheta,itime) + A4_imag*a3
!
!                    end if
!
!                end do !itime
!            end do !itheta
!        end do !iradius
!
!        !c1_real(:,2,2) = 1._rk
!
!
!        ! Compute 1-4 characteristics from extrapolation: difference in radius-local mean and boundary average
!        call self%compute_boundary_average(worker,bc_comm,density_bar,vel1_bar,vel2_bar,vel3_bar,pressure_bar, &
!                                                          density_avg,vel1_avg,vel2_avg,vel3_avg,pressure_avg)
!        c_avg = sqrt(gam*pressure_avg/density_avg)
!
!
!        ! Convert back to primitive variable perturbations
!        call characteristics_to_primitive(self,worker,bc_comm,          &
!                                          c1_real,       c1_imag,       &
!                                          c2_real,       c2_imag,       &
!                                          c3_real,       c3_imag,       &
!                                          c4_real,       c4_imag,       &
!                                          c5_real,       c5_imag,       &
!                                          density_real,  density_imag,  &
!                                          vel1_real,     vel1_imag,     &
!                                          vel2_real,     vel2_imag,     &
!                                          vel3_real,     vel3_imag,     &
!                                          pressure_real, pressure_imag)
!
!
!
!
!
!        ! Get boundary condition Total Temperature, Total Pressure, and normal vector
!        PT   = self%bcproperties%compute('Total Pressure',   worker%time(),worker%coords())
!        TT   = self%bcproperties%compute('Total Temperature',worker%time(),worker%coords())
!
!        ! Get user-input normal vector and normalize
!        n1 = self%bcproperties%compute('Normal-1', worker%time(), worker%coords())
!        n2 = self%bcproperties%compute('Normal-2', worker%time(), worker%coords())
!        n3 = self%bcproperties%compute('Normal-3', worker%time(), worker%coords())
!
!        !   Explicit allocation to handle GCC bug:
!        !       GCC/GFortran Bugzilla Bug 52162 
!        allocate(nmag(size(n1)), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!        nmag = sqrt(n1*n1 + n2*n2 + n3*n3)
!        n1 = n1/nmag
!        n2 = n2/nmag
!        n3 = n3/nmag
!
!
!        ! Override spatio-temporal mean according to specified total conditions
!        T_avg = TT(1)*(pressure_avg/PT(1))**((gam-ONE)/gam)
!        density_avg = pressure_avg/(T_avg*Rgas)
!        vmag = sqrt(TWO*cp*(TT(1)-T_avg))
!        vel1_avg = n1(1)*vmag
!        vel2_avg = n2(1)*vmag
!        vel3_avg = n3(1)*vmag
!
!        density_real(:,1,1)  = density_avg
!        vel1_real(:,1,1)     = vel1_avg
!        vel2_real(:,1,1)     = vel2_avg
!        vel3_real(:,1,1)     = vel3_avg
!        pressure_real(:,1,1) = pressure_avg
!
!
!
!    end subroutine compute_absorbing_inlet
!    !********************************************************************************


    !>  DANIEL'S FORMULATION
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_absorbing_inlet(self,worker,bc_comm, &
                                       density_real,    &
                                       density_imag,    &
                                       vel1_real,       &
                                       vel1_imag,       &
                                       vel2_real,       &
                                       vel2_imag,       &
                                       vel3_real,       &
                                       vel3_imag,       &
                                       pressure_real,   &
                                       pressure_imag)
        class(giles_HB_base_t),  intent(inout)   :: self
        type(chidg_worker_t),    intent(inout)   :: worker
        type(mpi_comm),          intent(in)      :: bc_comm
        type(AD_D),              intent(inout)   :: density_real(:,:,:)
        type(AD_D),              intent(inout)   :: density_imag(:,:,:)
        type(AD_D),              intent(inout)   :: vel1_real(:,:,:)
        type(AD_D),              intent(inout)   :: vel1_imag(:,:,:)
        type(AD_D),              intent(inout)   :: vel2_real(:,:,:)
        type(AD_D),              intent(inout)   :: vel2_imag(:,:,:)
        type(AD_D),              intent(inout)   :: vel3_real(:,:,:)
        type(AD_D),              intent(inout)   :: vel3_imag(:,:,:)
        type(AD_D),              intent(inout)   :: pressure_real(:,:,:)
        type(AD_D),              intent(inout)   :: pressure_imag(:,:,:)


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                                      &
            c1,    c2,    c3,    c4,    c5,                                             &
            c1_3d, c2_3d, c3_3d, c4_3d, c5_3d,                                          &
            c1_1d, c2_1d, c3_1d, c4_1d, c5_1d,                                          &
            c_bar, ddensity, dvel1, dvel2, dvel3, dpressure, expect_zero,               &
            density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar

        type(AD_D), allocatable, dimension(:,:,:) ::        &
            c1_real, c2_real, c3_real, c4_real, c5_real,    &
            c1_imag, c2_imag, c3_imag, c4_imag, c5_imag

        type(AD_D)  :: pressure_avg, vel1_avg, vel2_avg, vel3_avg, density_avg, c_avg, T_avg,       &
                       density_bar_r, vel1_bar_r, vel2_bar_r, vel3_bar_r, pressure_bar_r, c_bar_r,  &
                       beta, s_real, s_imag, s_arg, kz, k1, k2, k3, k4_real, k4_imag, k5_real, k5_imag, vmag,   &
                       a1_real, a2_real, a3_real, a4_real, a5_real, &
                       a1_imag, a2_imag, a3_imag, a4_imag, a5_imag, pyra

        real(rk),       allocatable, dimension(:)   :: PT, TT, n1, n2, n3, nmag, pitch
        real(rk)                                    :: theta_offset, omega, lm, a1, a2, a3, a4
        integer(ik)                                 :: iradius, igq, ierr, itheta, ntheta, itime, ntime

        pitch  = self%bcproperties%compute('Pitch',worker%time(),worker%coords())

        ! Get spatio-temporal average at radial stations
        density_bar  = density_real(:,1,1)
        vel1_bar     = vel1_real(:,1,1)
        vel2_bar     = vel2_real(:,1,1)
        vel3_bar     = vel3_real(:,1,1)
        pressure_bar = pressure_real(:,1,1)

        ntheta = size(c5_real,2)
        ntime  = size(c5_real,3)


!        ! Compute Fourier decomposition at set of radial stations: 
!        call self%primitive_to_characteristics(worker,bc_COMM,                  &
!                                               density_real,  density_imag,     &
!                                               vel1_real,     vel1_imag,        &
!                                               vel2_real,     vel2_imag,        &
!                                               vel3_real,     vel3_imag,        &
!                                               pressure_real, pressure_imag,    &
!                                               c1_real,       c1_imag,          &
!                                               c2_real,       c2_imag,          &
!                                               c3_real,       c3_imag,          &
!                                               c4_real,       c4_imag,          &
!                                               c5_real,       c5_imag)
!
!        ! Handle temporal average(steady) nonreflecting condition (:,:,1)
!        do iradius = 1,size(self%r)
!            ! Get average parts
!            density_bar_r  = density_bar(iradius)
!            vel1_bar_r     = vel1_bar(iradius)
!            vel2_bar_r     = vel2_bar(iradius)
!            vel3_bar_r     = vel3_bar(iradius)
!            pressure_bar_r = pressure_bar(iradius)
!            c_bar_r        = sqrt(gam*pressure_bar_r/density_bar_r)
!
!            ! starting with 2 here because the first mode is treated with 1D characteristics
!            do itheta = 2,ntheta
!                ! Account for sign(mode) in the calculation of beta. The second half of the
!                ! modes are negative frequencies.
!                if (itheta <= (ntheta-1)/2 + 1) then
!                    beta = sqrt(c_bar_r*c_bar_r  -  (vel3_bar_r*vel3_bar_r + vel2_bar_r*vel2_bar_r))
!                else if (itheta > (ntheta-1)/2 + 1) then
!                    beta = -sqrt(c_bar_r*c_bar_r  -  (vel3_bar_r*vel3_bar_r + vel2_bar_r*vel2_bar_r))
!                end if
!
!                ! The imaginary part of beta has already been accounted for in
!                ! the expressions for A2 and A3
!                A3_real = -TWO*vel3_bar_r*vel2_bar_r/(vel2_bar_r*vel2_bar_r + beta*beta)
!                A3_imag = -TWO*beta*vel3_bar_r/(vel2_bar_r*vel2_bar_r + beta*beta)
!
!                A4_real = (beta*beta - vel2_bar_r*vel2_bar_r)/(beta*beta + vel2_bar_r*vel2_bar_r)
!                A4_imag = -TWO*beta*vel2_bar_r/(beta*beta + vel2_bar_r*vel2_bar_r)
!
!                c1_real(iradius,itheta,1) = ZERO
!                c1_imag(iradius,itheta,1) = ZERO
!
!                c2_real(iradius,itheta,1) = ZERO
!                c2_imag(iradius,itheta,1) = ZERO
!
!                c3_real(iradius,itheta,1) = -(ONE/(c_bar_r + vel3_bar_r))*(vel2_bar_r*c5_real(iradius,itheta,1) - beta*c5_imag(iradius,itheta,1))
!                c3_imag(iradius,itheta,1) = -(ONE/(c_bar_r + vel3_bar_r))*(vel2_bar_r*c5_imag(iradius,itheta,1) + beta*c5_real(iradius,itheta,1))
!
!                c4_real(iradius,itheta,1) = (ONE/((c_bar_r+vel3_bar_r)**TWO))*((vel2_bar_r*vel2_bar_r - beta*beta)*c5_real(iradius,itheta,1) - TWO*beta*c5_imag(iradius,itheta,1))
!                c4_imag(iradius,itheta,1) = (ONE/((c_bar_r+vel3_bar_r)**TWO))*((vel2_bar_r*vel2_bar_r - beta*beta)*c5_imag(iradius,itheta,1) + TWO*beta*c5_real(iradius,itheta,1))
!
!            end do !itheta
!        end do !iradius


        ! Project
        print*, 'WARNING! HARDCODED ACCESS TO FIRST FREQUENCY ONLY!!!'
        print*, 'WARNING! CHECK DEFINITION OF lm PITCH!'
        do iradius = 1,size(self%r)
            ! Get average parts
            density_bar_r  = density_bar(iradius)
            vel1_bar_r     = vel1_bar(iradius)
            vel2_bar_r     = vel2_bar(iradius)
            vel3_bar_r     = vel3_bar(iradius)
            pressure_bar_r = pressure_bar(iradius)
            c_bar_r        = sqrt(gam*pressure_bar_r/density_bar_r)

            ! starting with 2 here because the first mode is treated with 1D characteristics
            do itheta = 1,ntheta
                do itime = 1,ntime

                    
                    ! Space-time average handled at the bottom
                    if (itime == 1 .and. itheta == 1) cycle

                    omega = worker%time_manager%freqs(1)
                    lm = TWO*PI*real(itheta-1,rk)/pitch(1)
                    kz = lm*c_bar_r/omega
                    pyra = (omega - kz*vel2_bar_r)**TWO - kz*kz*(c_bar_r**TWO - vel3_bar_r**TWO)

                    ! Compute k1,k2,k3
                    k1 = (omega - kz*vel2_bar_r)/vel3_bar_r
                    k2 = (omega - kz*vel2_bar_r)/vel3_bar_r
                    k3 = (omega - kz*vel2_bar_r)/vel3_bar_r

                    if (pyra >= ZERO) then
                        k4_real = (-vel3_bar_r*(omega - kz*vel2_bar_r) + c_bar_r*sqrt(pyra))/(c_bar_r**TWO - vel3_bar_r**TWO)
                        k4_imag = ZERO*k4_real
                        k5_real = (-vel3_bar_r*(omega - kz*vel2_bar_r) - c_bar_r*sqrt(pyra))/(c_bar_r**TWO - vel3_bar_r**TWO)
                        k5_imag = ZERO*k4_real
                    else
                        k4_real = (-vel3_bar_r*(omega - kz*vel2_bar_r))/(c_bar_r**TWO - vel3_bar_r**TWO)
                        k4_imag = c_bar_r*sqrt(pyra)/(c_bar_r**TWO - vel3_bar_r**TWO)
                        k5_real = (-vel3_bar_r*(omega - kz*vel2_bar_r))/(c_bar_r**TWO - vel3_bar_r**TWO)
                        k5_imag = -c_bar_r*sqrt(pyra)/(c_bar_r**TWO - vel3_bar_r**TWO)
                    end if


                    a1_real = (ONE/density_bar_r)*density_real(iradius,itheta,itime)  -  (ONE/(density_bar_r*c_bar_r*c_bar_r))*pressure_real(iradius,itheta,itime)
                    a1_imag = (ONE/density_bar_r)*density_imag(iradius,itheta,itime)  -  (ONE/(density_bar_r*c_bar_r*c_bar_r))*pressure_imag(iradius,itheta,itime)
                    
                    a2_real = (ONE/c_bar_r)*vel1_bar_r
                    a2_imag = (ONE/c_bar_r)*vel1_bar_r

                    a3_real = (-kz/(c_bar_r*(k1*k1 + kz*kz)))*vel3_real(iradius,itheta,itime)  +  (k1/(c_bar_r*(k1*k1+kz*kz)))*vel2_real(iradius,itheta,itime) - (kz/(density_bar_r*c_bar_r*vel3_bar_r*(k1*k1+kz*kz)))*pressure_real(iradius,itheta,itime)
                    a3_imag = (-kz/(c_bar_r*(k1*k1 + kz*kz)))*vel3_imag(iradius,itheta,itime)  +  (k1/(c_bar_r*(k1*k1+kz*kz)))*vel2_imag(iradius,itheta,itime) - (kz/(density_bar_r*c_bar_r*vel3_bar_r*(k1*k1+kz*kz)))*pressure_imag(iradius,itheta,itime)

                    ! Default, incoming waves set to zero
                    a1_real = ZERO
                    a1_imag = ZERO
                    a2_real = ZERO
                    a2_imag = ZERO
                    a3_real = ZERO
                    a3_imag = ZERO

                    ! User-specified amplitudes for incoming waves
                    if (itheta == 2 .and. itime == 2) then
                        a3_real = 1.
                    end if

                    ! Reset perturbation values
                    density_real(iradius,itheta,itime) = ZERO
                    density_imag(iradius,itheta,itime) = ZERO
                    vel1_real(iradius,itheta,itime) = ZERO
                    vel1_imag(iradius,itheta,itime) = ZERO
                    vel2_real(iradius,itheta,itime) = ZERO
                    vel2_imag(iradius,itheta,itime) = ZERO
                    vel3_real(iradius,itheta,itime) = ZERO
                    vel3_imag(iradius,itheta,itime) = ZERO
                    pressure_real(iradius,itheta,itime) = ZERO
                    pressure_imag(iradius,itheta,itime) = ZERO

                    ! Accumulate from absorbing amplitudes
                    vel3_real(iradius,itheta,itime) = -c_bar_r*kz*a3_real
                    vel2_real(iradius,itheta,itime) =  c_bar_r*k1*a3_real

                end do !itime
            end do !itheta
        end do !iradius


        ! Compute 1-4 characteristics from extrapolation: difference in radius-local mean and boundary average
        call self%compute_boundary_average(worker,bc_comm,density_bar,vel1_bar,vel2_bar,vel3_bar,pressure_bar, &
                                                          density_avg,vel1_avg,vel2_avg,vel3_avg,pressure_avg)
        c_avg = sqrt(gam*pressure_avg/density_avg)



        ! Get boundary condition Total Temperature, Total Pressure, and normal vector
        PT   = self%bcproperties%compute('Total Pressure',   worker%time(),worker%coords())
        TT   = self%bcproperties%compute('Total Temperature',worker%time(),worker%coords())

        ! Get user-input normal vector and normalize
        n1 = self%bcproperties%compute('Normal-1', worker%time(), worker%coords())
        n2 = self%bcproperties%compute('Normal-2', worker%time(), worker%coords())
        n3 = self%bcproperties%compute('Normal-3', worker%time(), worker%coords())

        !   Explicit allocation to handle GCC bug:
        !       GCC/GFortran Bugzilla Bug 52162 
        allocate(nmag(size(n1)), stat=ierr)
        if (ierr /= 0) call AllocationError

        nmag = sqrt(n1*n1 + n2*n2 + n3*n3)
        n1 = n1/nmag
        n2 = n2/nmag
        n3 = n3/nmag


        ! Override spatio-temporal mean according to specified total conditions
        T_avg = TT(1)*(pressure_avg/PT(1))**((gam-ONE)/gam)
        density_avg = pressure_avg/(T_avg*Rgas)
        vmag = sqrt(TWO*cp*(TT(1)-T_avg))
        vel1_avg = n1(1)*vmag
        vel2_avg = n2(1)*vmag
        vel3_avg = n3(1)*vmag

        density_real(:,1,1)  = density_avg
        vel1_real(:,1,1)     = vel1_avg
        vel2_real(:,1,1)     = vel2_avg
        vel3_real(:,1,1)     = vel3_avg
        pressure_real(:,1,1) = pressure_avg



    end subroutine compute_absorbing_inlet
    !********************************************************************************





    !>
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_absorbing_outlet(self,worker,bc_comm, &
                                        density_Fts_real,    &
                                        density_Fts_imag,    &
                                        vel1_Fts_real,       &
                                        vel1_Fts_imag,       &
                                        vel2_Fts_real,       &
                                        vel2_Fts_imag,       &
                                        vel3_Fts_real,       &
                                        vel3_Fts_imag,       &
                                        pressure_Fts_real,   &
                                        pressure_Fts_imag)
        class(giles_HB_base_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        type(mpi_comm),                             intent(in)      :: bc_comm
        type(AD_D),                                 intent(inout)   :: density_Fts_real(:,:,:)
        type(AD_D),                                 intent(inout)   :: density_Fts_imag(:,:,:)
        type(AD_D),                                 intent(inout)   :: vel1_Fts_real(:,:,:)
        type(AD_D),                                 intent(inout)   :: vel1_Fts_imag(:,:,:)
        type(AD_D),                                 intent(inout)   :: vel2_Fts_real(:,:,:)
        type(AD_D),                                 intent(inout)   :: vel2_Fts_imag(:,:,:)
        type(AD_D),                                 intent(inout)   :: vel3_Fts_real(:,:,:)
        type(AD_D),                                 intent(inout)   :: vel3_Fts_imag(:,:,:)
        type(AD_D),                                 intent(inout)   :: pressure_Fts_real(:,:,:)
        type(AD_D),                                 intent(inout)   :: pressure_Fts_imag(:,:,:)


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                                      &
            c1,    c2,    c3,    c4,    c5,                                             &
            c1_3d, c2_3d, c3_3d, c4_3d, c5_3d,                                          &
            c1_1d, c2_1d, c3_1d, c4_1d, c5_1d,                                          &
            c_bar, ddensity, dvel1, dvel2, dvel3, dpressure, expect_zero,               &
            density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar

        type(AD_D), allocatable, dimension(:,:,:) ::                            &
            c1_hat_real, c2_hat_real, c3_hat_real, c4_hat_real, c5_hat_real,    &
            c1_hat_imag, c2_hat_imag, c3_hat_imag, c4_hat_imag, c5_hat_imag

        type(AD_D)  :: pressure_avg, vel1_avg, vel2_avg, vel3_avg, density_avg, c_avg,              &
                       density_bar_r, vel1_bar_r, vel2_bar_r, vel3_bar_r, pressure_bar_r, c_bar_r,  &
                       A3_real, A3_imag, A4_real, A4_imag, beta, s_real, s_imag, s_arg, lambda

        type(point_t),  allocatable                 :: coords(:)
        real(rk),       allocatable, dimension(:)   :: p_user, r, pitch
        real(rk)                                    :: theta_offset, omega, lm
        integer(ik)                                 :: iradius, igq, ierr, itheta, ntheta, itime, ntime

        ! Get spatio-temporal average at radial stations
        density_bar  = density_Fts_real(:,1,1)
        vel1_bar     = vel1_Fts_real(:,1,1)
        vel2_bar     = vel2_Fts_real(:,1,1)
        vel3_bar     = vel3_Fts_real(:,1,1)
        pressure_bar = pressure_Fts_real(:,1,1)

        ! Retrieve target average pressure
        p_user = self%bcproperties%compute('Average Pressure',worker%time(),worker%coords())
        pitch  = self%bcproperties%compute('Pitch',worker%time(),worker%coords())


        ! Compute Fourier decomposition at set of radial stations: 
        call self%primitive_to_characteristics(worker,bc_COMM,                          &
                                               density_Fts_real,  density_Fts_imag,     &
                                               vel1_Fts_real,     vel1_Fts_imag,        &
                                               vel2_Fts_real,     vel2_Fts_imag,        &
                                               vel3_Fts_real,     vel3_Fts_imag,        &
                                               pressure_Fts_real, pressure_Fts_imag,    &
                                               c1_hat_real,       c1_hat_imag,          &
                                               c2_hat_real,       c2_hat_imag,          &
                                               c3_hat_real,       c3_hat_imag,          &
                                               c4_hat_real,       c4_hat_imag,          &
                                               c5_hat_real,       c5_hat_imag)

        ! Handle temporal average(steady) nonreflecting condition (:,:,1)
        ntheta = size(c5_hat_real,2)
        ntime  = size(c5_hat_real,3)
        do iradius = 1,size(self%r)
            ! Get average parts
            density_bar_r  = density_bar(iradius)
            vel1_bar_r     = vel1_bar(iradius)
            vel2_bar_r     = vel2_bar(iradius)
            vel3_bar_r     = vel3_bar(iradius)
            pressure_bar_r = pressure_bar(iradius)
            c_bar_r        = sqrt(gam*pressure_bar_r/density_bar_r)

            ! starting with 2 here because the first mode is treated with 1D characteristics
            do itheta = 2,ntheta
                ! Account for sign(mode) in the calculation of beta. The second half of the
                ! modes are negative frequencies.
                if (itheta <= (ntheta-1)/2 + 1) then
                    beta = sqrt(c_bar_r*c_bar_r  -  (vel3_bar_r*vel3_bar_r + vel2_bar_r*vel2_bar_r))
                else if (itheta > (ntheta-1)/2 + 1) then
                    beta = -sqrt(c_bar_r*c_bar_r  -  (vel3_bar_r*vel3_bar_r + vel2_bar_r*vel2_bar_r))
                end if

                ! The imaginary part of beta has already been accounted for in
                ! the expressions for A2 and A3
                A3_real = -TWO*vel3_bar_r*vel2_bar_r/(vel2_bar_r*vel2_bar_r + beta*beta)
                A3_imag = -TWO*beta*vel3_bar_r/(vel2_bar_r*vel2_bar_r + beta*beta)

                A4_real = (beta*beta - vel2_bar_r*vel2_bar_r)/(beta*beta + vel2_bar_r*vel2_bar_r)
                A4_imag = -TWO*beta*vel2_bar_r/(beta*beta + vel2_bar_r*vel2_bar_r)

                c5_hat_real(iradius,itheta,1) = (A3_real*c3_hat_real(iradius,itheta,1) - A3_imag*c3_hat_imag(iradius,itheta,1))  &   ! A3*c3 (real)
                                              - (A4_real*c4_hat_real(iradius,itheta,1) - A4_imag*c4_hat_imag(iradius,itheta,1))      ! A4*c4 (real)
                c5_hat_imag(iradius,itheta,1) = (A3_imag*c3_hat_real(iradius,itheta,1) + A3_real*c3_hat_imag(iradius,itheta,1))  &   ! A3*c3 (imag)
                                              - (A4_imag*c4_hat_real(iradius,itheta,1) + A4_real*c4_hat_imag(iradius,itheta,1))      ! A4*c4 (imag)
            end do !itheta
        end do !iradius


        ! Handle unsteady nonreflecting condition (:,:,2:)
        print*, 'WARNING! HARDCODED ACCESS TO FIRST FREQUENCY ONLY!!!'
        print*, 'WARNING! CHECK DEFINITION OF lm PITCH!'
        do iradius = 1,size(self%r)
            ! Get average parts
            density_bar_r  = density_bar(iradius)
            vel1_bar_r     = vel1_bar(iradius)
            vel2_bar_r     = vel2_bar(iradius)
            vel3_bar_r     = vel3_bar(iradius)
            pressure_bar_r = pressure_bar(iradius)
            c_bar_r        = sqrt(gam*pressure_bar_r/density_bar_r)


            ! starting with 2 here because the first mode is treated with 1D characteristics
            do itheta = 1,ntheta
                do itime = 2,ntime

                    omega = worker%time_manager%freqs(1)
                    lm = TWO*PI*real(itheta-1,rk)/pitch(1)
                    lambda = lm*c_bar_r/omega
                    s_arg = ONE - (c_bar_r*c_bar_r - vel3_bar_r*vel3_bar_r)*lambda*lambda/((c_bar_r - vel2_bar_r*lambda)**TWO)

                    if (s_arg >= ZERO) then
                        s_real = sqrt(s_arg)
                        s_imag = ZERO*s_arg
                    else
                        s_real = ZERO*s_arg
                        s_imag = sqrt(-s_arg) ! minus because we have already factored out sqrt(-1) = i into the formula below
                    end if

                    A3_real = TWO*vel3_bar_r*lambda*(ONE + s_real)/((c_bar_r-vel2_bar_r*lambda)*((ONE+s_real)**TWO + s_imag*s_imag))
                    A3_imag = -TWO*vel3_bar_r*lambda*s_imag/((c_bar_r-vel2_bar_r*lambda)*((ONE+s_real)**TWO + s_imag*s_imag))

                    A4_real = (ONE - s_real*s_real - s_imag*s_imag)/((ONE+s_real)**TWO + s_imag*s_imag)
                    A4_imag = -TWO*s_imag/((ONE+s_real)**TWO + s_imag*s_imag)

                    c5_hat_real(iradius,itheta,itime) = (A3_real*c3_hat_real(iradius,itheta,itime) - A3_imag*c3_hat_imag(iradius,itheta,itime))  &   ! A3*c3 (real)
                                                      + (A4_real*c4_hat_real(iradius,itheta,itime) - A4_imag*c4_hat_imag(iradius,itheta,itime))      ! A4*c4 (real)
                    c5_hat_imag(iradius,itheta,itime) = (A3_imag*c3_hat_real(iradius,itheta,itime) + A3_real*c3_hat_imag(iradius,itheta,itime))  &   ! A3*c3 (imag)
                                                      + (A4_imag*c4_hat_real(iradius,itheta,itime) + A4_real*c4_hat_imag(iradius,itheta,itime))      ! A4*c4 (imag)

                end do !itime
            end do !itheta
        end do !iradius


        ! Compute 1-4 characteristics from extrapolation: difference in radius-local mean and boundary average
        call self%compute_boundary_average(worker,bc_comm,density_bar,vel1_bar,vel2_bar,vel3_bar,pressure_bar, &
                                                          density_avg,vel1_avg,vel2_avg,vel3_avg,pressure_avg)
        c_avg = sqrt(gam*pressure_avg/density_avg)


        ! Compute contribution due to driving of spatio-temporal average
        ddensity  = density_bar  - density_avg 
        dvel1     = vel1_bar     - vel1_avg
        dvel2     = vel2_bar     - vel2_avg
        dvel3     = vel3_bar     - vel3_avg
        dpressure = pressure_bar - pressure_avg
        c1_1d = ZERO*ddensity
        c2_1d = ZERO*ddensity
        c3_1d = ZERO*ddensity
        c4_1d = ZERO*ddensity
        c5_1d = ZERO*ddensity
        do iradius = 1,size(self%r)
            c1_1d(iradius) = -c_avg*c_avg*ddensity(iradius)    +  dpressure(iradius)
            c2_1d(iradius) = density_avg*c_avg*dvel1(iradius)
            c3_1d(iradius) = density_avg*c_avg*dvel2(iradius)
            c4_1d(iradius) = density_avg*c_avg*dvel3(iradius)  +  dpressure(iradius)
            c5_1d(iradius) = -TWO*(pressure_avg - p_user(1))
        end do
        ! Drive the spatio-temporal average
        c1_hat_real(:,1,1) = c1_hat_real(:,1,1) + c1_1d(:)
        c2_hat_real(:,1,1) = c2_hat_real(:,1,1) + c2_1d(:)
        c3_hat_real(:,1,1) = c3_hat_real(:,1,1) + c3_1d(:)
        c4_hat_real(:,1,1) = c4_hat_real(:,1,1) + c4_1d(:)
        c5_hat_real(:,1,1) = c5_hat_real(:,1,1) + c5_1d(:)



        call characteristics_to_primitive(self,worker,bc_comm,                  &
                                          c1_hat_real,       c1_hat_imag,       &
                                          c2_hat_real,       c2_hat_imag,       &
                                          c3_hat_real,       c3_hat_imag,       &
                                          c4_hat_real,       c4_hat_imag,       &
                                          c5_hat_real,       c5_hat_imag,       &
                                          density_Fts_real,  density_Fts_imag,  &
                                          vel1_Fts_real,     vel1_Fts_imag,     &
                                          vel2_Fts_real,     vel2_Fts_imag,     &
                                          vel3_Fts_real,     vel3_Fts_imag,     &
                                          pressure_Fts_real, pressure_Fts_imag)


    end subroutine compute_absorbing_outlet
    !********************************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine primitive_to_characteristics(self,worker,bc_comm,                    &
                                            density_Fts_real,  density_Fts_imag,    &
                                            vel1_Fts_real,     vel1_Fts_imag,       &
                                            vel2_Fts_real,     vel2_Fts_imag,       &
                                            vel3_Fts_real,     vel3_Fts_imag,       &
                                            pressure_Fts_real, pressure_Fts_imag,   &
                                            c1_hat_real,       c1_hat_imag,         &
                                            c2_hat_real,       c2_hat_imag,         &
                                            c3_hat_real,       c3_hat_imag,         &
                                            c4_hat_real,       c4_hat_imag,         &
                                            c5_hat_real,       c5_hat_imag)
        class(giles_HB_base_t),         intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        type(mpi_comm),                 intent(in)      :: bc_comm
        type(AD_D),                     intent(in)      :: density_Fts_real(:,:,:)
        type(AD_D),                     intent(in)      :: density_Fts_imag(:,:,:)
        type(AD_D),                     intent(in)      :: vel1_Fts_real(:,:,:)
        type(AD_D),                     intent(in)      :: vel1_Fts_imag(:,:,:)
        type(AD_D),                     intent(in)      :: vel2_Fts_real(:,:,:)
        type(AD_D),                     intent(in)      :: vel2_Fts_imag(:,:,:)
        type(AD_D),                     intent(in)      :: vel3_Fts_real(:,:,:)
        type(AD_D),                     intent(in)      :: vel3_Fts_imag(:,:,:)
        type(AD_D),                     intent(in)      :: pressure_Fts_real(:,:,:)
        type(AD_D),                     intent(in)      :: pressure_Fts_imag(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c1_hat_real(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c1_hat_imag(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c2_hat_real(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c2_hat_imag(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c3_hat_real(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c3_hat_imag(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c4_hat_real(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c4_hat_imag(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c5_hat_real(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c5_hat_imag(:,:,:)

        type(AD_D), allocatable,    dimension(:)    ::                                          &
            density, mom1, mom2, mom3, energy, vel1, vel2, vel3, pressure,                      &
            c1_real_tmp,      c2_real_tmp,   c3_real_tmp,   c4_real_tmp,   c5_real_tmp,         &
            c1_imag_tmp,      c2_imag_tmp,   c3_imag_tmp,   c4_imag_tmp,   c5_imag_tmp,         &
            density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar

        type(AD_D)  :: c_bar

        integer(ik)             :: nmodes, nradius, iradius, itheta, ntheta, imode, ierr, ntime, itime
        real(rk)                :: shift_r, shift_i
        real(rk),   allocatable :: pitch(:)

        ! Get spatio-temporal average at radial stations
        density_bar  = density_Fts_real(:,1,1)
        vel1_bar     = vel1_Fts_real(:,1,1)
        vel2_bar     = vel2_Fts_real(:,1,1)
        vel3_bar     = vel3_Fts_real(:,1,1)
        pressure_bar = pressure_Fts_real(:,1,1)

        ! Define Fourier discretization
        nmodes  = self%nfourier_space
        ntheta  = 1 + (nmodes-1)*2
        nradius = size(self%r)
        ntime = worker%time_manager%ntime

        pitch  = self%bcproperties%compute('Pitch',time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])

        ! Allocate storage in result
        allocate(c1_hat_real(nradius,ntheta,ntime), c1_hat_imag(nradius,ntheta,ntime),  &
                 c2_hat_real(nradius,ntheta,ntime), c2_hat_imag(nradius,ntheta,ntime),  &
                 c3_hat_real(nradius,ntheta,ntime), c3_hat_imag(nradius,ntheta,ntime),  &
                 c4_hat_real(nradius,ntheta,ntime), c4_hat_imag(nradius,ntheta,ntime),  &
                 c5_hat_real(nradius,ntheta,ntime), c5_hat_imag(nradius,ntheta,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError
        c1_hat_real = ZERO*density_bar(1)
        c2_hat_real = ZERO*density_bar(1)
        c3_hat_real = ZERO*density_bar(1)
        c4_hat_real = ZERO*density_bar(1)
        c5_hat_real = ZERO*density_bar(1)
        c1_hat_imag = ZERO*density_bar(1)
        c2_hat_imag = ZERO*density_bar(1)
        c3_hat_imag = ZERO*density_bar(1)
        c4_hat_imag = ZERO*density_bar(1)
        c5_hat_imag = ZERO*density_bar(1)

        ! Convert Fourier modes of primitive varibles to 1D characteristics
        do iradius = 1,nradius
            c_bar = sqrt(gam*pressure_bar(iradius)/density_bar(iradius))
            do itheta = 1,ntheta
                do itime = 1,ntime
                    c1_hat_real(iradius,itheta,itime) = -(c_bar*c_bar)*density_Fts_real(iradius,itheta,itime)             +  (ONE)*pressure_Fts_real(iradius,itheta,itime)
                    c2_hat_real(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel1_Fts_real(iradius,itheta,itime)
                    c3_hat_real(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel2_Fts_real(iradius,itheta,itime)
                    c4_hat_real(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel3_Fts_real(iradius,itheta,itime)  +  (ONE)*pressure_Fts_real(iradius,itheta,itime)
                    c5_hat_real(iradius,itheta,itime) = -(density_bar(iradius)*c_bar)*vel3_Fts_real(iradius,itheta,itime) +  (ONE)*pressure_Fts_real(iradius,itheta,itime)
                                               
                    c1_hat_imag(iradius,itheta,itime) = -(c_bar*c_bar)*density_Fts_imag(iradius,itheta,itime)             +  (ONE)*pressure_Fts_imag(iradius,itheta,itime)
                    c2_hat_imag(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel1_Fts_imag(iradius,itheta,itime)
                    c3_hat_imag(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel2_Fts_imag(iradius,itheta,itime)
                    c4_hat_imag(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel3_Fts_imag(iradius,itheta,itime)  +  (ONE)*pressure_Fts_imag(iradius,itheta,itime)
                    c5_hat_imag(iradius,itheta,itime) = -(density_bar(iradius)*c_bar)*vel3_Fts_imag(iradius,itheta,itime) +  (ONE)*pressure_Fts_imag(iradius,itheta,itime)
                end do !itime
            end do !itheta
        end do !iradius

    end subroutine primitive_to_characteristics
    !*********************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine characteristics_to_primitive(self,worker,bc_comm,                  &
                                            c1_hat_real,       c1_hat_imag,       &
                                            c2_hat_real,       c2_hat_imag,       &
                                            c3_hat_real,       c3_hat_imag,       &
                                            c4_hat_real,       c4_hat_imag,       &
                                            c5_hat_real,       c5_hat_imag,       &
                                            density_Fts_real,  density_Fts_imag,  &
                                            vel1_Fts_real,     vel1_Fts_imag,     &
                                            vel2_Fts_real,     vel2_Fts_imag,     &
                                            vel3_Fts_real,     vel3_Fts_imag,     &
                                            pressure_Fts_real, pressure_Fts_imag)
        class(giles_HB_base_t),         intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        type(mpi_comm),                 intent(in)      :: bc_comm
        type(AD_D),                     intent(in)      :: c1_hat_real(:,:,:)
        type(AD_D),                     intent(in)      :: c1_hat_imag(:,:,:)
        type(AD_D),                     intent(in)      :: c2_hat_real(:,:,:)
        type(AD_D),                     intent(in)      :: c2_hat_imag(:,:,:)
        type(AD_D),                     intent(in)      :: c3_hat_real(:,:,:)
        type(AD_D),                     intent(in)      :: c3_hat_imag(:,:,:)
        type(AD_D),                     intent(in)      :: c4_hat_real(:,:,:)
        type(AD_D),                     intent(in)      :: c4_hat_imag(:,:,:)
        type(AD_D),                     intent(in)      :: c5_hat_real(:,:,:)
        type(AD_D),                     intent(in)      :: c5_hat_imag(:,:,:)
        type(AD_D),                     intent(inout)   :: density_Fts_real(:,:,:)
        type(AD_D),                     intent(inout)   :: density_Fts_imag(:,:,:)
        type(AD_D),                     intent(inout)   :: vel1_Fts_real(:,:,:)
        type(AD_D),                     intent(inout)   :: vel1_Fts_imag(:,:,:)
        type(AD_D),                     intent(inout)   :: vel2_Fts_real(:,:,:)
        type(AD_D),                     intent(inout)   :: vel2_Fts_imag(:,:,:)
        type(AD_D),                     intent(inout)   :: vel3_Fts_real(:,:,:)
        type(AD_D),                     intent(inout)   :: vel3_Fts_imag(:,:,:)
        type(AD_D),                     intent(inout)   :: pressure_Fts_real(:,:,:)
        type(AD_D),                     intent(inout)   :: pressure_Fts_imag(:,:,:)

        type(AD_D), allocatable,    dimension(:)    ::  &
            density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar
            
        type(AD_D)  :: density_bar_r, pressure_bar_r, c_bar_r

        integer(ik)             :: iradius, itheta, itime

        ! Get spatio-temporal average at radial stations
        density_bar  = density_Fts_real(:,1,1)
        vel1_bar     = vel1_Fts_real(:,1,1)
        vel2_bar     = vel2_Fts_real(:,1,1)
        vel3_bar     = vel3_Fts_real(:,1,1)
        pressure_bar = pressure_Fts_real(:,1,1)

        ! Convert characteristic Fourier modes back to primitive Fourier modes and store
        do iradius = 1,size(c1_hat_real,1)
            ! Get radius-local averages
            density_bar_r  = density_bar(iradius)
            pressure_bar_r = pressure_bar(iradius)
            c_bar_r        = sqrt(gam*pressure_bar_r/density_bar_r)
            do itheta = 1,size(c1_hat_real,2)
                do itime = 1,size(c1_hat_real,3)
                    density_Fts_real(iradius,itheta,itime) = -(ONE/(c_bar_r*c_bar_r))*c1_hat_real(iradius,itheta,itime) + (ONE/(TWO*c_bar_r*c_bar_r))*c4_hat_real(iradius,itheta,itime) + (ONE/(TWO*c_bar_r*c_bar_r))*c5_hat_real(iradius,itheta,itime)
                    vel1_Fts_real(iradius,itheta,itime) = (ONE/(density_bar_r*c_bar_r))*c2_hat_real(iradius,itheta,itime)
                    vel2_Fts_real(iradius,itheta,itime) = (ONE/(density_bar_r*c_bar_r))*c3_hat_real(iradius,itheta,itime)
                    vel3_Fts_real(iradius,itheta,itime) = (ONE/(TWO*density_bar_r*c_bar_r))*c4_hat_real(iradius,itheta,itime) - (ONE/(TWO*density_bar_r*c_bar_r))*c5_hat_real(iradius,itheta,itime)
                    pressure_Fts_real(iradius,itheta,itime) = HALF*c4_hat_real(iradius,itheta,itime) + HALF*c5_hat_real(iradius,itheta,itime)

                    density_Fts_imag(iradius,itheta,itime) = -(ONE/(c_bar_r*c_bar_r))*c1_hat_imag(iradius,itheta,itime) + (ONE/(TWO*c_bar_r*c_bar_r))*c4_hat_imag(iradius,itheta,itime) + (ONE/(TWO*c_bar_r*c_bar_r))*c5_hat_imag(iradius,itheta,itime)
                    vel1_Fts_imag(iradius,itheta,itime) = (ONE/(density_bar_r*c_bar_r))*c2_hat_imag(iradius,itheta,itime)
                    vel2_Fts_imag(iradius,itheta,itime) = (ONE/(density_bar_r*c_bar_r))*c3_hat_imag(iradius,itheta,itime)
                    vel3_Fts_imag(iradius,itheta,itime) = (ONE/(TWO*density_bar_r*c_bar_r))*c4_hat_imag(iradius,itheta,itime) - (ONE/(TWO*density_bar_r*c_bar_r))*c5_hat_imag(iradius,itheta,itime)
                    pressure_Fts_imag(iradius,itheta,itime) = HALF*c4_hat_imag(iradius,itheta,itime) + HALF*c5_hat_imag(iradius,itheta,itime)
                end do
            end do
        end do 

    end subroutine characteristics_to_primitive
    !*********************************************************************************











    !> Compute boundary average by averaging spatio-temporal average over
    !! radial stations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/16/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_boundary_average(self,worker,bc_comm,density_bar,vel1_bar,vel2_bar,vel3_bar,pressure_bar, &
                                                            density_avg,vel1_avg,vel2_avg,vel3_avg,pressure_avg)
        class(giles_HB_base_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(in)      :: worker
        type(mpi_comm),                             intent(in)      :: bc_comm
        type(AD_D),                                 intent(in)      :: density_bar(:)
        type(AD_D),                                 intent(in)      :: vel1_bar(:)
        type(AD_D),                                 intent(in)      :: vel2_bar(:)
        type(AD_D),                                 intent(in)      :: vel3_bar(:)
        type(AD_D),                                 intent(in)      :: pressure_bar(:)
        type(AD_D),                                 intent(inout)   :: density_avg
        type(AD_D),                                 intent(inout)   :: vel1_avg
        type(AD_D),                                 intent(inout)   :: vel2_avg
        type(AD_D),                                 intent(inout)   :: vel3_avg
        type(AD_D),                                 intent(inout)   :: pressure_avg

        real(rk)    :: area, dr
        integer(ik) :: irad


        density_avg  = ZERO*density_bar(1)
        vel1_avg     = ZERO*density_bar(1)
        vel2_avg     = ZERO*density_bar(1)
        vel3_avg     = ZERO*density_bar(1)
        pressure_avg = ZERO*density_bar(1)

        if (worker%coordinate_system() == 'Cartesian') then
            dr = self%r(2) - self%r(1)
            area = ZERO
            do irad = 1,size(density_bar)-1
                density_avg  = density_avg  + dr*(density_bar(irad+1) +density_bar(irad))/TWO
                vel1_avg     = vel1_avg     + dr*(vel1_bar(irad+1)    +vel1_bar(irad))/TWO
                vel2_avg     = vel2_avg     + dr*(vel2_bar(irad+1)    +vel2_bar(irad))/TWO
                vel3_avg     = vel3_avg     + dr*(vel3_bar(irad+1)    +vel3_bar(irad))/TWO
                pressure_avg = pressure_avg + dr*(pressure_bar(irad+1)+pressure_bar(irad))/TWO
                area = area + dr
            end do

        else if (worker%coordinate_system() == 'Cylindrical') then
            dr = self%r(2) - self%r(1)
            area = ZERO
            do irad = 1,size(density_bar)-1
                density_avg  = density_avg  + dr*(self%r(irad+1)*density_bar(irad+1)  + self%r(irad)*density_bar(irad))/TWO
                vel1_avg     = vel1_avg     + dr*(self%r(irad+1)*vel1_bar(irad+1)     + self%r(irad)*vel1_bar(irad))/TWO
                vel2_avg     = vel2_avg     + dr*(self%r(irad+1)*vel2_bar(irad+1)     + self%r(irad)*vel2_bar(irad))/TWO
                vel3_avg     = vel3_avg     + dr*(self%r(irad+1)*vel3_bar(irad+1)     + self%r(irad)*vel3_bar(irad))/TWO
                pressure_avg = pressure_avg + dr*(self%r(irad+1)*pressure_bar(irad+1) + self%r(irad)*pressure_bar(irad))/TWO
                area = area + dr*(self%r(irad+1)+self%r(irad))/TWO
            end do

        end if

        density_avg  = density_avg  / area
        vel1_avg     = vel1_avg     / area
        vel2_avg     = vel2_avg     / area
        vel3_avg     = vel3_avg     / area
        pressure_avg = pressure_avg / area

    end subroutine compute_boundary_average
    !********************************************************************************





    !>  Determine rmin, rmax, and average theta.
    !!
    !!  Initialize radial stations and reference theta for Fourier transform.
    !!  
    !!  Defines:
    !!      self%r(:)
    !!      self%theta_ref
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/26/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine analyze_bc_geometry(self,mesh,group_ID,bc_comm)
        class(giles_HB_base_t),  intent(inout)   :: self
        type(mesh_t),                               intent(in)      :: mesh
        integer(ik),                                intent(in)      :: group_ID
        type(mpi_comm),                             intent(in)      :: bc_comm

        integer(ik) :: bc_idomain_l, bc_ielement_l, patch_ID, face_ID, iface, ierr, inode
        real(rk)    :: face_rmin, face_rmax, local_rmin, local_rmax, global_rmin, global_rmax,  &
                       face_thetamin, face_thetamax, local_thetamin, local_thetamax, global_thetamin, global_thetamax
        real(rk)    :: ref_nodes(8,3), physical_nodes(8,3)

        ! Search for min/max radius on local processor
        local_rmin     =  HUGE(1._rk) ! any radius will be smaller than this, so it is guarunteed to be reset.
        local_rmax     = -HUGE(1._rk) ! any radius will be larger than this, so it is guarunteed to be reset.
        local_thetamin =  HUGE(1._rk) ! any theta will be smaller than this, so it is guarunteed to be reset.
        local_thetamax = -HUGE(1._rk) ! any theta will be larger than this, so it is guarunteed to be reset.
        do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
            do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                iface = mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)

                ! Pick points to evaluate coordinates
                if (iface == XI_MIN) then
                    ref_nodes(1,:) = [-ONE, -ONE, -ONE]
                    ref_nodes(2,:) = [-ONE, -ONE, ZERO]
                    ref_nodes(3,:) = [-ONE, -ONE,  ONE]
                    ref_nodes(4,:) = [-ONE, ZERO,  ONE]
                    ref_nodes(5,:) = [-ONE,  ONE,  ONE]
                    ref_nodes(6,:) = [-ONE,  ONE, ZERO]
                    ref_nodes(7,:) = [-ONE,  ONE, -ONE]
                    ref_nodes(8,:) = [-ONE, ZERO, -ONE]
                else if (iface == XI_MAX) then
                    ref_nodes(1,:) = [ONE, -ONE, -ONE]
                    ref_nodes(2,:) = [ONE, -ONE, ZERO]
                    ref_nodes(3,:) = [ONE, -ONE,  ONE]
                    ref_nodes(4,:) = [ONE, ZERO,  ONE]
                    ref_nodes(5,:) = [ONE,  ONE,  ONE]
                    ref_nodes(6,:) = [ONE,  ONE, ZERO]
                    ref_nodes(7,:) = [ONE,  ONE, -ONE]
                    ref_nodes(8,:) = [ONE, ZERO, -ONE]

                else if (iface == ETA_MIN) then
                    ref_nodes(1,:) = [ -ONE, -ONE, -ONE]
                    ref_nodes(2,:) = [ -ONE, -ONE, ZERO]
                    ref_nodes(3,:) = [ -ONE, -ONE,  ONE]
                    ref_nodes(4,:) = [ ZERO, -ONE,  ONE]
                    ref_nodes(5,:) = [  ONE, -ONE,  ONE]
                    ref_nodes(6,:) = [  ONE, -ONE, ZERO]
                    ref_nodes(7,:) = [  ONE, -ONE, -ONE]
                    ref_nodes(8,:) = [ ZERO, -ONE, -ONE]

                else if (iface == ETA_MAX) then
                    ref_nodes(1,:) = [ -ONE, ONE, -ONE]
                    ref_nodes(2,:) = [ -ONE, ONE, ZERO]
                    ref_nodes(3,:) = [ -ONE, ONE,  ONE]
                    ref_nodes(4,:) = [ ZERO, ONE,  ONE]
                    ref_nodes(5,:) = [  ONE, ONE,  ONE]
                    ref_nodes(6,:) = [  ONE, ONE, ZERO]
                    ref_nodes(7,:) = [  ONE, ONE, -ONE]
                    ref_nodes(8,:) = [ ZERO, ONE, -ONE]

                else if (iface == ZETA_MIN) then
                    ref_nodes(1,:) = [ -ONE, -ONE, -ONE]
                    ref_nodes(2,:) = [ -ONE, ZERO, -ONE]
                    ref_nodes(3,:) = [ -ONE,  ONE, -ONE]
                    ref_nodes(4,:) = [ ZERO,  ONE, -ONE]
                    ref_nodes(5,:) = [  ONE,  ONE, -ONE]
                    ref_nodes(6,:) = [  ONE, ZERO, -ONE]
                    ref_nodes(7,:) = [  ONE, -ONE, -ONE]
                    ref_nodes(8,:) = [ ZERO, -ONE, -ONE]

                else if (iface == ZETA_MAX) then
                    ref_nodes(1,:) = [ -ONE, -ONE, ONE]
                    ref_nodes(2,:) = [ -ONE, ZERO, ONE]
                    ref_nodes(3,:) = [ -ONE,  ONE, ONE]
                    ref_nodes(4,:) = [ ZERO,  ONE, ONE]
                    ref_nodes(5,:) = [  ONE,  ONE, ONE]
                    ref_nodes(6,:) = [  ONE, ZERO, ONE]
                    ref_nodes(7,:) = [  ONE, -ONE, ONE]
                    ref_nodes(8,:) = [ ZERO, -ONE, ONE]

                else
                    call chidg_signal(FATAL,"outlet_giles_quasi3d_unsteady_HB: analyze_bc_geometry, invalid face indec.")
                end if

                ! Evaluate physical coordinates on face edges
                bc_idomain_l  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
                bc_ielement_l = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                do inode = 1,size(ref_nodes,1)
                    physical_nodes(inode,:) = mesh%domain(bc_idomain_l)%elems(bc_ielement_l)%physical_point(ref_nodes(inode,:),'Deformed')
                end do !inode

                ! Get face min/max radius
                face_rmin     = minval(physical_nodes(:,1))
                face_rmax     = maxval(physical_nodes(:,1))
                face_thetamin = minval(physical_nodes(:,2))
                face_thetamax = maxval(physical_nodes(:,2))

                ! Update processor-local value if new min/max values were found on the face
                if (face_rmin < local_rmin)         local_rmin     = face_rmin
                if (face_rmax > local_rmax)         local_rmax     = face_rmax
                if (face_thetamin < local_thetamin) local_thetamin = face_thetamin
                if (face_thetamax > local_thetamax) local_thetamax = face_thetamax

            end do !face_ID
        end do !patch_ID

        ! Reduce processor local values to determine boundary-global min/max values
        call MPI_AllReduce(local_rmin,    global_rmin,    1,MPI_REAL8,MPI_MIN,bc_comm,ierr)
        call MPI_AllReduce(local_rmax,    global_rmax,    1,MPI_REAL8,MPI_MAX,bc_comm,ierr)
        call MPI_AllReduce(local_thetamin,global_thetamin,1,MPI_REAL8,MPI_MIN,bc_comm,ierr)
        call MPI_AllReduce(local_thetamax,global_thetamax,1,MPI_REAL8,MPI_MAX,bc_comm,ierr)
        
        ! Create radial stations
        self%r = linspace(global_rmin,global_rmax,self%nr)

        ! Compute theta_ref
        self%theta_ref = (global_thetamin + global_thetamax)/TWO

    end subroutine analyze_bc_geometry
    !********************************************************************************




    !>  For each radial station, a theta-discretization exists upon which the 
    !!  discrete Fourier transform operation is computed. Since the theta-discretization
    !!  spans the entire boundary, a general interpolation procedure is used to fill
    !!  solution values that could come from many different elements on the boundary.
    !!
    !!  This procedure finds donor elements for each node in the theta-discretization
    !!  at each radial station. Additionally, the location of the physical coordinate
    !!  within the donor element's reference coordinate system is also found.
    !!
    !!  This procedure
    !!
    !!  Initializes:
    !!  ------------------------------
    !!  self%theta(:)
    !!  self%donor(:,:)
    !!  self%donor_node(:,:,:)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/28/2018
    !!
    !----------------------------------------------------------------------------------
    subroutine initialize_fourier_discretization(self,mesh,group_ID,bc_comm)
        class(giles_HB_base_t),  intent(inout)   :: self
        type(mesh_t),                               intent(in)      :: mesh
        integer(ik),                                intent(in)      :: group_ID
        type(mpi_comm),                             intent(in)      :: bc_comm

        integer(ik)                 :: nmodes, ncoeff, nradius, ntheta, idomain_l, ielement_l, iface, &
                                       iradius, itheta, ierr, noverset
        real(rk)                    :: dtheta, dtheta_n, midpoint(3), try_offset(3), node(3), z
        real(rk),       allocatable :: pitch(:)
        character(:),   allocatable :: user_msg
        logical                     :: donor_found


        ! Determine z-location of some face on the boundary and assume 
        ! entire boundary is constant-z
        idomain_l  = mesh%bc_patch_group(group_ID)%patch(1)%idomain_l()
        ielement_l = mesh%bc_patch_group(group_ID)%patch(1)%ielement_l(1)
        iface      = mesh%bc_patch_group(group_ID)%patch(1)%iface(1)
        if (iface == XI_MIN) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([-ONE,ZERO,ZERO],'Deformed')
        else if (iface == XI_MAX) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ ONE,ZERO,ZERO],'Deformed')
        else if (iface == ETA_MIN) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ZERO,-ONE,ZERO],'Deformed')
        else if (iface == ETA_MAX) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ZERO, ONE,ZERO],'Deformed')
        else if (iface == ZETA_MIN) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ZERO,ZERO,-ONE],'Deformed')
        else if (iface == ZETA_MAX) then
            midpoint = mesh%domain(idomain_l)%elems(ielement_l)%physical_point([ZERO,ZERO, ONE],'Deformed')
        end if
        z = midpoint(3)

        ! Define Fourier discretization
        nmodes  = self%nfourier_space
        ncoeff  = 1 + (nmodes-1)*2
        nradius = size(self%r)
        ntheta  = ncoeff
        
        ! Initialize theta discretization parameters
        pitch  = self%bcproperties%compute('Pitch',time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])
        dtheta = pitch(1)
        dtheta_n = dtheta/ntheta

        ! Construct theta discretization at each radius
        allocate(self%theta(size(self%r),ntheta), stat=ierr)
        if (ierr /= 0) call AllocationError
        do itheta = 1,ntheta
            self%theta(:,itheta) = self%theta_ref + (itheta-1)*dtheta_n
        end do

        ! Donor search offset, if needed
        try_offset = [ZERO, -pitch(1), ZERO]

        ! For each radial station, initialized donor for each node in theta grid
        allocate(self%donor(nradius,ntheta), self%donor_node(nradius,ntheta,3), stat=ierr)
        if (ierr /= 0) call AllocationError
        do iradius = 1,size(self%r)
            noverset = 0
            do itheta = 1,ntheta

                node = [self%r(iradius), self%theta(iradius,itheta), z]

                ! Try processor-LOCAL elements
                call find_gq_donor(mesh,                                    &
                                   node,                                    &
                                   [ZERO,ZERO,ZERO],                        &
                                   face_info_constructor(0,0,0,0,0),        &   ! we don't really have a receiver face
                                   self%donor(iradius,itheta),              &
                                   self%donor_node(iradius,itheta,1:3),     &
                                   donor_found)

                ! Try LOCAL elements with try_offset if still not found 
                if ( .not. donor_found ) then
                    call find_gq_donor(mesh,                                    &
                                       node,                                    &
                                       try_offset,                              &
                                       face_info_constructor(0,0,0,0,0),        &   ! we don't really have a receiver face
                                       self%donor(iradius,itheta),              &
                                       self%donor_node(iradius,itheta,1:3),     &
                                       donor_found)
                    if (donor_found) then
                        noverset=noverset+1
                        self%theta(iradius,itheta) = self%theta(iradius,itheta) - pitch(1)
                    end if
                end if

                ! Try PARALLEL_ELEMENTS if donor not found amongst local elements
                if (.not. donor_found) then
                    call find_gq_donor_parallel(mesh,                                   &
                                                node,                                   &
                                                [ZERO,ZERO,ZERO],                       &
                                                face_info_constructor(0,0,0,0,0),       &   ! we don't really have a receiver face
                                                self%donor(iradius,itheta),             &
                                                self%donor_node(iradius,itheta,1:3),    &
                                                donor_found)
                end if

                
                ! Try PARALLEL_ELEMENTS with try_offset if still not found 
                if ( .not. donor_found ) then
                    call find_gq_donor_parallel(mesh,                                   &
                                                node,                                   &
                                                try_offset,                             &
                                                face_info_constructor(0,0,0,0,0),       &   ! we don't really have a receiver face
                                                self%donor(iradius,itheta),             &
                                                self%donor_node(iradius,itheta,1:3),    &
                                                donor_found)
                    if (donor_found) then
                        noverset=noverset+1
                        self%theta(iradius,itheta) = self%theta(iradius,itheta) - pitch(1)
                    end if
                end if 


                ! Abort if we didn't find a donor
                user_msg = "bc_giles_HB_base%initialize_fourier_discretization: &
                            no donor element found for Fourier discretization node."
                if (.not. donor_found) call chidg_signal(FATAL,user_msg)

            end do !itheta

            ! Shift arrays so that we start with the theta_min point
            self%donor(iradius,:)        = cshift(self%donor(iradius,:), -noverset, dim=1)
            self%donor_node(iradius,:,:) = cshift(self%donor_node(iradius,:,:), -noverset, dim=1)
            self%theta(iradius,:)        = cshift(self%theta(iradius,:), -noverset, dim=1)

        end do !iradius

    end subroutine initialize_fourier_discretization
    !**************************************************************************************




end module bc_giles_HB_base
