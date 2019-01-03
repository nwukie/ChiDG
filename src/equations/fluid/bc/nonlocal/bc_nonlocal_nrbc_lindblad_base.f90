module bc_nonlocal_nrbc_lindblad_base
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, FOUR, HALF, ME, CYLINDRICAL,    &
                                      XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, PI
    use mod_io,                 only: verbosity
    use mod_fluid,              only: Rgas, cp, gam
    use mod_interpolation,      only: interpolate_linear, interpolate_linear_ad
    use mod_gridspace,          only: linspace
    use mod_dft,                only: dft, idft_eval
    use mod_chimera,            only: find_gq_donor, find_gq_donor_parallel, multi_donor_rule_t

    use type_point,             only: point_t
    use type_mesh,              only: mesh_t
    use type_bc_state,          only: bc_state_t
    use type_bc_patch,          only: bc_patch_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_face_info,         only: face_info_t, face_info_constructor
    use type_element_info,      only: element_info_t
    use type_ivector,           only: ivector_t
    use mod_chidg_mpi,          only: IRANK
    use mod_interpolate,        only: interpolate_face_autodiff
    use mpi_f08,                only: MPI_REAL8, MPI_INTEGER4, MPI_AllReduce, mpi_comm, MPI_BCast, MPI_MIN, MPI_MAX, MPI_SUM
    use ieee_arithmetic,        only: ieee_is_nan
    use DNAD_D
    implicit none



    !>  Base class for Fourier interface for turbomachinery blade rows. Can be used
    !!  for boundary conditions or as an interface.
    !!
    !!  ****
    !!  ASSUMES CONSTANT Z-PLANE 
    !!  ****
    !!
    !!  References:
    !!              
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2018
    !!
    !---------------------------------------------------------------------------------
    type, public, abstract, extends(bc_state_t) :: nonlocal_nrbc_lindblad_base_t

        integer(ik) :: nr = 10
        integer(ik) :: nfourier_space = 20

        real(rk),               allocatable :: r(:)
        real(rk)                            :: theta_ref = ZERO

        integer(ik)                         :: nfaces_a         ! Global count across all procs in group
        real(rk),               allocatable :: theta_a(:,:)     ! (nr,ntheta)
        real(rk)                            :: theta_ref_a
        real(rk)                            :: z_ref_a
        type(element_info_t),   allocatable :: donor_a(:,:)
        real(rk),               allocatable :: donor_node_a(:,:,:)

        integer(ik)                         :: nfaces_b         ! Global count across all procs in group
        real(rk),               allocatable :: theta_b(:,:)     ! (nr,ntheta)
        real(rk)                            :: theta_ref_b
        real(rk)                            :: z_ref_b
        type(element_info_t),   allocatable :: donor_b(:,:)
        real(rk),               allocatable :: donor_node_b(:,:,:)

    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: init_bc_coupling     ! Initialize global coupling
        procedure   :: init_bc_postcomm     ! Implement specialized initialization

        procedure   :: get_q_side

        procedure   :: compute_temporal_dft
        procedure   :: compute_spatial_dft
        procedure   :: compute_spatial_idft_gq
        procedure   :: compute_temporal_idft_gq
        procedure   :: interpolate_raux_to_rgq
        procedure   :: primitive_to_eigenmodes
        procedure   :: eigenmodes_to_primitive
        procedure   :: compute_eigenvalues
        procedure   :: compute_boundary_average
        procedure   :: analyze_bc_geometry
        procedure   :: initialize_fourier_discretization

        procedure   :: get_face_side

        procedure   :: primitive_to_characteristics
        procedure   :: characteristics_to_primitive

    end type nonlocal_nrbc_lindblad_base_t
    !*********************************************************************************


    type, extends(multi_donor_rule_t), public :: rule_prefers_A
    contains
        procedure, nopass :: select_donor => select_donor_A
    end type

    type, extends(multi_donor_rule_t), public :: rule_prefers_B
    contains
        procedure, nopass :: select_donor => select_donor_B
    end type

contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(nonlocal_nrbc_lindblad_base_t),   intent(inout) :: self

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
        class(nonlocal_nrbc_lindblad_base_t), intent(inout)   :: self
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
        class(nonlocal_nrbc_lindblad_base_t),     intent(inout)   :: self
        type(mesh_t),               intent(inout)   :: mesh
        integer(ik),                intent(in)      :: group_ID
        type(mpi_comm),             intent(in)      :: bc_comm

        call self%analyze_bc_geometry(mesh,group_ID,bc_comm,side='A')
        call self%analyze_bc_geometry(mesh,group_ID,bc_comm,side='B')

        call self%initialize_fourier_discretization(mesh,group_ID,bc_comm,side='A')
        call self%initialize_fourier_discretization(mesh,group_ID,bc_comm,side='B')

    end subroutine init_bc_postcomm
    !*************************************************************************************



    
    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/25/2018
    !!
    !------------------------------------------------------------------------------------
    subroutine get_q_side(self,worker,bc_comm,side,  &
                          density,                   &
                          vel1,                      &
                          vel2,                      &
                          vel3,                      &
                          pressure)
        class(nonlocal_nrbc_lindblad_base_t),         intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        type(mpi_comm),                 intent(in)      :: bc_comm
        character(1),                   intent(in)      :: side
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
        if (side=='A') ntheta  = size(self%theta_a,2)
        if (side=='B') ntheta  = size(self%theta_b,2)
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
                if (side=='A') then
                    density(iradius,:,itime) = worker%interpolate_field_general('Density',    donors=self%donor_a(iradius,:), donor_nodes=self%donor_node_a(iradius,:,:), itime=itime)
                    mom1(iradius,:,itime)    = worker%interpolate_field_general('Momentum-1', donors=self%donor_a(iradius,:), donor_nodes=self%donor_node_a(iradius,:,:), itime=itime)
                    mom2(iradius,:,itime)    = worker%interpolate_field_general('Momentum-2', donors=self%donor_a(iradius,:), donor_nodes=self%donor_node_a(iradius,:,:), itime=itime)
                    mom3(iradius,:,itime)    = worker%interpolate_field_general('Momentum-3', donors=self%donor_a(iradius,:), donor_nodes=self%donor_node_a(iradius,:,:), itime=itime)
                    energy(iradius,:,itime)  = worker%interpolate_field_general('Energy',     donors=self%donor_a(iradius,:), donor_nodes=self%donor_node_a(iradius,:,:), itime=itime)
                else if (side=='B') then
                    density(iradius,:,itime) = worker%interpolate_field_general('Density',    donors=self%donor_b(iradius,:), donor_nodes=self%donor_node_b(iradius,:,:), itime=itime)
                    mom1(iradius,:,itime)    = worker%interpolate_field_general('Momentum-1', donors=self%donor_b(iradius,:), donor_nodes=self%donor_node_b(iradius,:,:), itime=itime)
                    mom2(iradius,:,itime)    = worker%interpolate_field_general('Momentum-2', donors=self%donor_b(iradius,:), donor_nodes=self%donor_node_b(iradius,:,:), itime=itime)
                    mom3(iradius,:,itime)    = worker%interpolate_field_general('Momentum-3', donors=self%donor_b(iradius,:), donor_nodes=self%donor_node_b(iradius,:,:), itime=itime)
                    energy(iradius,:,itime)  = worker%interpolate_field_general('Energy',     donors=self%donor_b(iradius,:), donor_nodes=self%donor_node_b(iradius,:,:), itime=itime)
                end if

                if (worker%coordinate_system() == 'Cylindrical') then
                    mom2(iradius,:,itime) = mom2(iradius,:,itime)/self%r(iradius)  ! convert to tangential momentum
                end if
            end do
        end do

        ! Compute velocities and pressure at each time
        vel1 = mom1/density
        vel2 = mom2/density
        vel3 = mom3/density
        pressure = (gam-ONE)*(energy - HALF*(mom1*mom1 + mom2*mom2 + mom3*mom3)/density)


    end subroutine get_q_side
    !************************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/12/2018
    !!
    !------------------------------------------------------------------------------------
    subroutine compute_temporal_dft(self,worker,bc_comm,                    &
                                    density_grid, vel1_grid, vel2_grid, vel3_grid, pressure_grid, c_grid, &
                                    density_Ft_real,  density_Ft_imag,      &
                                    vel1_Ft_real,     vel1_Ft_imag,         &
                                    vel2_Ft_real,     vel2_Ft_imag,         &
                                    vel3_Ft_real,     vel3_Ft_imag,         &
                                    pressure_Ft_real, pressure_Ft_imag,     &
                                    c_Ft_real,        c_Ft_imag)
        class(nonlocal_nrbc_lindblad_base_t),         intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        type(mpi_comm),                 intent(in)      :: bc_comm
        type(AD_D),     allocatable,    intent(inout)   :: density_grid(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: vel1_grid(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: vel2_grid(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: vel3_grid(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: pressure_grid(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c_grid(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: density_Ft_real(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: density_Ft_imag(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: vel1_Ft_real(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: vel1_Ft_imag(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: vel2_Ft_real(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: vel2_Ft_imag(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: vel3_Ft_real(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: vel3_Ft_imag(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: pressure_Ft_real(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: pressure_Ft_imag(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c_Ft_real(:,:,:)
        type(AD_D),     allocatable,    intent(inout)   :: c_Ft_imag(:,:,:)


        type(AD_D), allocatable,    dimension(:)    ::                                          &
            density_real_tmp, vel1_real_tmp, vel2_real_tmp, vel3_real_tmp, pressure_real_tmp, c_real_tmp,   &
            density_imag_tmp, vel1_imag_tmp, vel2_imag_tmp, vel3_imag_tmp, pressure_imag_tmp, c_imag_tmp

        type(AD_D)  :: density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, c_bar

        integer(ik) :: nradius, ntheta, iradius, itheta, imode, itime, ntime, ierr


        ! Define Fourier space discretization to determine
        ! number of theta-samples are being taken
        !ntheta  = size(self%theta,2)
        !nradius = size(self%r)
        nradius = size(density_grid,1)
        ntheta  = size(density_grid,2)
        ntime   = worker%mesh%domain(worker%element_info%idomain_l)%elems(worker%element_info%ielement_l)%ntime

        ! Allocate storage for temporal dft
        allocate(density_Ft_real( nradius,ntheta,ntime), density_Ft_imag( nradius,ntheta,ntime),  &
                 vel1_Ft_real(    nradius,ntheta,ntime), vel1_Ft_imag(    nradius,ntheta,ntime),  &
                 vel2_Ft_real(    nradius,ntheta,ntime), vel2_Ft_imag(    nradius,ntheta,ntime),  &
                 vel3_Ft_real(    nradius,ntheta,ntime), vel3_Ft_imag(    nradius,ntheta,ntime),  &
                 pressure_Ft_real(nradius,ntheta,ntime), pressure_Ft_imag(nradius,ntheta,ntime),  &
                 c_Ft_real(       nradius,ntheta,ntime), c_Ft_imag(       nradius,ntheta,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Perform Fourier decomposition at each radial station.
        do iradius = 1,nradius
            do itheta = 1,ntheta
                call dft(density_grid(iradius,itheta,:),  ZERO*density_grid(iradius,itheta,:),  density_real_tmp,  density_imag_tmp )
                call dft(vel1_grid(iradius,itheta,:),     ZERO*vel1_grid(iradius,itheta,:),     vel1_real_tmp,     vel1_imag_tmp    )
                call dft(vel2_grid(iradius,itheta,:),     ZERO*vel2_grid(iradius,itheta,:),     vel2_real_tmp,     vel2_imag_tmp    )
                call dft(vel3_grid(iradius,itheta,:),     ZERO*vel3_grid(iradius,itheta,:),     vel3_real_tmp,     vel3_imag_tmp    )
                call dft(pressure_grid(iradius,itheta,:), ZERO*pressure_grid(iradius,itheta,:), pressure_real_tmp, pressure_imag_tmp)
                call dft(c_grid(iradius,itheta,:),        ZERO*c_grid(iradius,itheta,:),        c_real_tmp,        c_imag_tmp       )

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
                c_Ft_real(       iradius,itheta,:) = c_real_tmp
                c_Ft_imag(       iradius,itheta,:) = c_imag_tmp

            end do !itheta
        end do !iradius

    end subroutine compute_temporal_dft
    !********************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_spatial_dft(self,worker,bc_comm,side,              &
                                   density_Ft_real,  density_Ft_imag,     &
                                   vel1_Ft_real,     vel1_Ft_imag,        &
                                   vel2_Ft_real,     vel2_Ft_imag,        &
                                   vel3_Ft_real,     vel3_Ft_imag,        &
                                   pressure_Ft_real, pressure_Ft_imag,    &
                                   c_Ft_real,        c_Ft_imag,           &
                                   density_Fts_real,  density_Fts_imag,   &
                                   vel1_Fts_real,     vel1_Fts_imag,      &
                                   vel2_Fts_real,     vel2_Fts_imag,      &
                                   vel3_Fts_real,     vel3_Fts_imag,      &
                                   pressure_Fts_real, pressure_Fts_imag,  &
                                   c_Fts_real,        c_Fts_imag)
        class(nonlocal_nrbc_lindblad_base_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        type(mpi_comm),                             intent(in)      :: bc_comm
        character(*),                               intent(in)      :: side
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
        type(AD_D),     allocatable,                intent(inout)   :: c_Ft_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c_Ft_imag(:,:,:)
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
        type(AD_D),     allocatable,                intent(inout)   :: c_Fts_real(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c_Fts_imag(:,:,:)

        type(AD_D), allocatable,    dimension(:)    ::                                          &
            density_real_tmp, vel1_real_tmp, vel2_real_tmp, vel3_real_tmp, pressure_real_tmp, c_real_tmp, &
            density_imag_tmp, vel1_imag_tmp, vel2_imag_tmp, vel3_imag_tmp, pressure_imag_tmp, c_imag_tmp

        integer(ik) :: nradius, ntheta, iradius, m, imode, itime, ntime, ierr
        real(rk)    :: shift_r, shift_i, shift_sign
        logical     :: negate_dft
        real(rk),   allocatable :: spatial_periodicity(:)

        ! Define Fourier space discretization to determine
        ! number of theta-samples being taken
        nradius = size(density_Ft_real,1)
        ntheta = size(density_Ft_real,2)
        ntime   = worker%mesh%domain(worker%element_info%idomain_l)%elems(worker%element_info%ielement_l)%ntime
        !spatial_periodicity = self%bcproperties%compute('Spatial Periodicity', time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])
        if (side=='A') then
            spatial_periodicity = self%bcproperties%compute('Pitch A', time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])
        else if (side=='B') then
            spatial_periodicity = self%bcproperties%compute('Pitch B', time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])
        end if
        
        ! Allocate storage in result
        allocate(density_Fts_real( nradius,ntheta,ntime), density_Fts_imag( nradius,ntheta,ntime),  &
                 vel1_Fts_real(    nradius,ntheta,ntime), vel1_Fts_imag(    nradius,ntheta,ntime),  &
                 vel2_Fts_real(    nradius,ntheta,ntime), vel2_Fts_imag(    nradius,ntheta,ntime),  &
                 vel3_Fts_real(    nradius,ntheta,ntime), vel3_Fts_imag(    nradius,ntheta,ntime),  &
                 pressure_Fts_real(nradius,ntheta,ntime), pressure_Fts_imag(nradius,ntheta,ntime),  &
                 c_Fts_real(       nradius,ntheta,ntime), c_Fts_imag(       nradius,ntheta,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Perform Fourier decomposition at each radial station.
        do iradius = 1,nradius
            do itime = 1,ntime

                ! this is consistens with Lindblad's formulation, which is used for the unsteady modes(itime>1)
                negate_dft = .true. 
                shift_sign = -ONE


                ! DFT in space
                call dft(density_Ft_real( iradius,:,itime), density_Ft_imag( iradius,:,itime), density_real_tmp,  density_imag_tmp,  negate=negate_dft)
                call dft(vel1_Ft_real(    iradius,:,itime), vel1_Ft_imag(    iradius,:,itime), vel1_real_tmp,     vel1_imag_tmp,     negate=negate_dft)
                call dft(vel2_Ft_real(    iradius,:,itime), vel2_Ft_imag(    iradius,:,itime), vel2_real_tmp,     vel2_imag_tmp,     negate=negate_dft)
                call dft(vel3_Ft_real(    iradius,:,itime), vel3_Ft_imag(    iradius,:,itime), vel3_real_tmp,     vel3_imag_tmp,     negate=negate_dft)
                call dft(pressure_Ft_real(iradius,:,itime), pressure_Ft_imag(iradius,:,itime), pressure_real_tmp, pressure_imag_tmp, negate=negate_dft)
                call dft(c_Ft_real(       iradius,:,itime), c_Ft_imag(       iradius,:,itime), c_real_tmp,        c_imag_tmp,        negate=negate_dft)


                ! Adjust Fourier coefficients so their phase is relative to self%theta_ref
                ! instead of the minimum theta of the transform.
                !
                !       q(relative to theta_ref) = q(relative to theta_min) * e^(j 2pi imode delta_theta/spatial_periodicity)
                !
                ! NOTE: self%theta(:,1) are defined to be the DFT-theta_min at each radius
                !
                do imode = 1,size(density_real_tmp)
                    m = get_lm(imode,size(density_real_tmp))
                    if (side=='A') then
                        shift_r = realpart(exp(shift_sign*cmplx(ZERO,ONE)*real(m,rk)*TWO*PI*(self%theta_ref-self%theta_a(iradius,1))/spatial_periodicity(1)))
                        shift_i = imagpart(exp(shift_sign*cmplx(ZERO,ONE)*real(m,rk)*TWO*PI*(self%theta_ref-self%theta_a(iradius,1))/spatial_periodicity(1)))
                    else if (side=='B') then
                        shift_r = realpart(exp(shift_sign*cmplx(ZERO,ONE)*real(m,rk)*TWO*PI*(self%theta_ref-self%theta_b(iradius,1))/spatial_periodicity(1)))
                        shift_i = imagpart(exp(shift_sign*cmplx(ZERO,ONE)*real(m,rk)*TWO*PI*(self%theta_ref-self%theta_b(iradius,1))/spatial_periodicity(1)))
                    else
                        call chidg_signal(FATAL,"giles_HB_base_t%compute_spatial_dft: Invalid input for argument 'side'. 'A' or 'B'.")
                    end if

                    density_Fts_real( iradius,imode,itime) = density_real_tmp(imode)*shift_r  - density_imag_tmp(imode)*shift_i
                    vel1_Fts_real(    iradius,imode,itime) = vel1_real_tmp(imode)*shift_r     - vel1_imag_tmp(imode)*shift_i
                    vel2_Fts_real(    iradius,imode,itime) = vel2_real_tmp(imode)*shift_r     - vel2_imag_tmp(imode)*shift_i
                    vel3_Fts_real(    iradius,imode,itime) = vel3_real_tmp(imode)*shift_r     - vel3_imag_tmp(imode)*shift_i
                    pressure_Fts_real(iradius,imode,itime) = pressure_real_tmp(imode)*shift_r - pressure_imag_tmp(imode)*shift_i
                    c_Fts_real(       iradius,imode,itime) = c_real_tmp(imode)*shift_r        - c_imag_tmp(imode)*shift_i

                    density_Fts_imag( iradius,imode,itime) = density_imag_tmp(imode)*shift_r  + density_real_tmp(imode)*shift_i
                    vel1_Fts_imag(    iradius,imode,itime) = vel1_imag_tmp(imode)*shift_r     + vel1_real_tmp(imode)*shift_i
                    vel2_Fts_imag(    iradius,imode,itime) = vel2_imag_tmp(imode)*shift_r     + vel2_real_tmp(imode)*shift_i
                    vel3_Fts_imag(    iradius,imode,itime) = vel3_imag_tmp(imode)*shift_r     + vel3_real_tmp(imode)*shift_i
                    pressure_Fts_imag(iradius,imode,itime) = pressure_imag_tmp(imode)*shift_r + pressure_real_tmp(imode)*shift_i
                    c_Fts_imag(       iradius,imode,itime) = c_imag_tmp(imode)*shift_r        + c_real_tmp(imode)*shift_i

                end do !imode

            end do !itime
        end do !iradius

        call write_line('WARNING: need scaling since dft is only over a single passage.', io_proc=IRANK, silence=(verbosity<5))
        call write_line('WARNING: check correct pitch in phase shift.', io_proc=IRANK, silence=(verbosity<5))

    end subroutine compute_spatial_dft
    !*********************************************************************************



    !>  Compute inverse Fourier transform of spatio-temporal Fourier coefficients
    !!  to quadrature nodes for the current face.
    !!
    !!  Input: q_hat(r_gq)
    !!  Output: q_check(r_gq,theta_gq)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/18/2018
    !!
    !----------------------------------------------------------------------------
    subroutine compute_spatial_idft_gq(self,worker,bc_comm,side, &
                                       density_hat_real,         &
                                       density_hat_imag,         &
                                       vel1_hat_real,            &
                                       vel1_hat_imag,            &
                                       vel2_hat_real,            &
                                       vel2_hat_imag,            &
                                       vel3_hat_real,            &
                                       vel3_hat_imag,            &
                                       pressure_hat_real,        &
                                       pressure_hat_imag,        &
                                       density_check_real,       &
                                       density_check_imag,       &
                                       vel1_check_real,          &
                                       vel1_check_imag,          &
                                       vel2_check_real,          &
                                       vel2_check_imag,          &
                                       vel3_check_real,          &
                                       vel3_check_imag,          &
                                       pressure_check_real,      &
                                       pressure_check_imag)
        class(nonlocal_nrbc_lindblad_base_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(mpi_comm),             intent(in)      :: bc_comm
        character(1),               intent(in)      :: side
        type(AD_D), allocatable,    intent(inout)   :: density_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: density_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: density_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: density_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_check_imag(:,:)

        type(point_t),  allocatable, dimension(:)   :: coords
        real(rk),       allocatable, dimension(:)   :: pitch
        real(rk)                                    :: theta_offset
        integer(ik) :: igq, itime
        logical :: negate_dft

        ! Get BC properties
        if (side == 'A') then
            pitch = self%bcproperties%compute('Pitch A', worker%time(),worker%coords())
        else if (side == 'B') then
            pitch = self%bcproperties%compute('Pitch B', worker%time(),worker%coords())
        end if

        ! Get gq coordinates 
        coords = worker%coords()

        ! Inverse DFT of spatio-temporal Fourier modes to give temporal Fourier
        ! modes at quadrature nodes.
        density_check_real  = ZERO*density_hat_real(:,1,:)
        vel1_check_real     = ZERO*density_hat_real(:,1,:)
        vel2_check_real     = ZERO*density_hat_real(:,1,:)
        vel3_check_real     = ZERO*density_hat_real(:,1,:)
        pressure_check_real = ZERO*density_hat_real(:,1,:)
        density_check_imag  = ZERO*density_hat_real(:,1,:)
        vel1_check_imag     = ZERO*density_hat_real(:,1,:)
        vel2_check_imag     = ZERO*density_hat_real(:,1,:)
        vel3_check_imag     = ZERO*density_hat_real(:,1,:)
        pressure_check_imag = ZERO*density_hat_real(:,1,:)

        
        negate_dft = .true.     ! this is consistent with Lindblad's formulation
        do igq = 1,size(coords)
            do itime = 1,size(density_hat_real,3)

                theta_offset = coords(igq)%c2_ - self%theta_ref
                ! **** WARNING: probably want ipdft_eval here ****
                call idft_eval(density_hat_real(igq,:,itime),       &
                               density_hat_imag(igq,:,itime),       &
                               [theta_offset]/pitch(1),             &
                               density_check_real(igq:igq,itime),   &
                               density_check_imag(igq:igq,itime),negate=negate_dft)

                call idft_eval(vel1_hat_real(igq,:,itime),          &
                               vel1_hat_imag(igq,:,itime),          &
                               [theta_offset]/pitch(1),             &
                               vel1_check_real(igq:igq,itime),      &
                               vel1_check_imag(igq:igq,itime),negate=negate_dft)

                call idft_eval(vel2_hat_real(igq,:,itime),          &
                               vel2_hat_imag(igq,:,itime),          &
                               [theta_offset]/pitch(1),             &
                               vel2_check_real(igq:igq,itime),      &
                               vel2_check_imag(igq:igq,itime),negate=negate_dft)

                call idft_eval(vel3_hat_real(igq,:,itime),          &
                               vel3_hat_imag(igq,:,itime),          &
                               [theta_offset]/pitch(1),             &
                               vel3_check_real(igq:igq,itime),      &
                               vel3_check_imag(igq:igq,itime),negate=negate_dft)

                call idft_eval(pressure_hat_real(igq,:,itime),      &
                               pressure_hat_imag(igq,:,itime),      &
                               [theta_offset]/pitch(1),             &
                               pressure_check_real(igq:igq,itime),  &
                               pressure_check_imag(igq:igq,itime),negate=negate_dft)
            end do !itime
        end do !igq

    end subroutine compute_spatial_idft_gq
    !***************************************************************************








    !>  Compute inverse Fourier transform of spatio-temporal Fourier coefficients
    !!  to quadrature nodes for the current face.
    !!
    !!  Input: q_check(r_gq,theta_gq)
    !!  Output: q(r_gq,theta_gq,itime)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/18/2018
    !!
    !----------------------------------------------------------------------------
    subroutine compute_temporal_idft_gq(self,worker,bc_comm,    &
                                        density_check_real,     &
                                        density_check_imag,     &
                                        vel1_check_real,        &
                                        vel1_check_imag,        &
                                        vel2_check_real,        &
                                        vel2_check_imag,        &
                                        vel3_check_real,        &
                                        vel3_check_imag,        &
                                        pressure_check_real,    &
                                        pressure_check_imag,    &
                                        density, vel1, vel2, vel3, pressure)
        class(nonlocal_nrbc_lindblad_base_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(mpi_comm),             intent(in)      :: bc_comm
        type(AD_D), allocatable,    intent(inout)   :: density_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: density_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_check_real(:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_check_imag(:,:)
        type(AD_D), allocatable,    intent(inout)   :: density(:)
        type(AD_D), allocatable,    intent(inout)   :: vel1(:)
        type(AD_D), allocatable,    intent(inout)   :: vel2(:)
        type(AD_D), allocatable,    intent(inout)   :: vel3(:)
        type(AD_D), allocatable,    intent(inout)   :: pressure(:)

        type(AD_D),     allocatable, dimension(:)   :: expect_zero, density_tmp, vel1_tmp, vel2_tmp, vel3_tmp, pressure_tmp
        integer(ik) :: igq

        ! Allocate storage for boundary state
        expect_zero = [AD_D(1)]
        density  = ZERO*density_check_real(:,1)
        vel1     = ZERO*density_check_real(:,1)
        vel2     = ZERO*density_check_real(:,1)
        vel3     = ZERO*density_check_real(:,1)
        pressure = ZERO*density_check_real(:,1)

        ! Inverse DFT of temporal Fourier modes to give primitive variables
        ! at quarature nodes for the current time instance.
        density_tmp  = [ZERO*density_check_real(1,1)]
        vel1_tmp     = [ZERO*density_check_real(1,1)]
        vel2_tmp     = [ZERO*density_check_real(1,1)]
        vel3_tmp     = [ZERO*density_check_real(1,1)]
        pressure_tmp = [ZERO*density_check_real(1,1)]
        do igq = 1,size(density_check_real,1)
            ! **** WARNING: probably want ipdft_eval here ****
            call idft_eval(density_check_real(igq,:),                                       &
                           density_check_imag(igq,:),                                       &
                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
                           density_tmp,                                                     &
                           expect_zero,symmetric=.true.)
            if (abs(expect_zero(1)) > 0.0000001) print*, 'WARNING: inverse transform returning complex values.'

            call idft_eval(vel1_check_real(igq,:),                                          &
                           vel1_check_imag(igq,:),                                          &
                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
                           vel1_tmp,                                                        &
                           expect_zero,symmetric=.true.)
            if (abs(expect_zero(1)) > 0.0000001) print*, 'WARNING: inverse transform returning complex values.'

            call idft_eval(vel2_check_real(igq,:),                                          &
                           vel2_check_imag(igq,:),                                          &
                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
                           vel2_tmp,                                                        &
                           expect_zero,symmetric=.true.)
            if (abs(expect_zero(1)) > 0.0000001) print*, 'WARNING: inverse transform returning complex values.'

            call idft_eval(vel3_check_real(igq,:),                                          &
                           vel3_check_imag(igq,:),                                          &
                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
                           vel3_tmp,                                                        &
                           expect_zero,symmetric=.true.)
            if (abs(expect_zero(1)) > 0.0000001) print*, 'WARNING: inverse transform returning complex values.'

            call idft_eval(pressure_check_real(igq,:),                                      &
                           pressure_check_imag(igq,:),                                      &
                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
                           pressure_tmp,                                                    &
                           expect_zero,symmetric=.true.)
            if (abs(expect_zero(1)) > 0.0000001) print*, 'WARNING: inverse transform returning complex values.'

            ! Accumulate total contribution from unsteady modes
            density(igq)  = density(igq)  + density_tmp(1)
            vel1(igq)     = vel1(igq)     + vel1_tmp(1)
            vel2(igq)     = vel2(igq)     + vel2_tmp(1)
            vel3(igq)     = vel3(igq)     + vel3_tmp(1)
            pressure(igq) = pressure(igq) + pressure_tmp(1)

        end do

    end subroutine compute_temporal_idft_gq
    !****************************************************************************








    !>  Interpolate Fourier coefficients from auxiliary grid to quadrature grid.
    !!
    !!  q_hat(r_gq) = I(q_hat(r_aux))
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/18/2018
    !!
    !-------------------------------------------------------------------------------------
    subroutine interpolate_raux_to_rgq(self,worker,bc_comm,                         &
                                       density_hat_real,     density_hat_imag,      &
                                       vel1_hat_real,        vel1_hat_imag,         &
                                       vel2_hat_real,        vel2_hat_imag,         &
                                       vel3_hat_real,        vel3_hat_imag,         &
                                       pressure_hat_real,    pressure_hat_imag,     &
                                       density_hat_real_gq,  density_hat_imag_gq,   &
                                       vel1_hat_real_gq,     vel1_hat_imag_gq,      &
                                       vel2_hat_real_gq,     vel2_hat_imag_gq,      &
                                       vel3_hat_real_gq,     vel3_hat_imag_gq,      &
                                       pressure_hat_real_gq, pressure_hat_imag_gq)
        class(nonlocal_nrbc_lindblad_base_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(mpi_comm),             intent(in)      :: bc_comm
        type(AD_D), allocatable,    intent(in)      :: density_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: density_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: vel1_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: vel1_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: vel2_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: vel2_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: vel3_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: vel3_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: pressure_hat_real(:,:,:)
        type(AD_D), allocatable,    intent(in)      :: pressure_hat_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: density_hat_real_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: density_hat_imag_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_hat_real_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_hat_imag_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_hat_real_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_hat_imag_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_hat_real_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_hat_imag_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_hat_real_gq(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_hat_imag_gq(:,:,:)
        
        type(point_t),  allocatable, dimension(:)   :: coords
        integer(ik) :: igq, itheta, itime, ierr


        ! Interpolate spatio-temporal Fourier coefficients to quadrature nodes
        ! linear interpolation between radial coordinates.
        coords = worker%coords() ! get gq coordinates for current face
        allocate(density_hat_real_gq( size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 density_hat_imag_gq( size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 vel1_hat_real_gq(    size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 vel1_hat_imag_gq(    size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 vel2_hat_real_gq(    size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 vel2_hat_imag_gq(    size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 vel3_hat_real_gq(    size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 vel3_hat_imag_gq(    size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 pressure_hat_real_gq(size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 pressure_hat_imag_gq(size(coords),size(density_hat_real,2),size(density_hat_real,3)), &
                 stat=ierr)
        if (ierr /= 0) call AllocationError
        density_hat_real_gq(:,:,:)  = ZERO*density_hat_real(1,1,1)
        density_hat_imag_gq(:,:,:)  = ZERO*density_hat_real(1,1,1)
        vel1_hat_real_gq(:,:,:)     = ZERO*density_hat_real(1,1,1)
        vel1_hat_imag_gq(:,:,:)     = ZERO*density_hat_real(1,1,1)
        vel2_hat_real_gq(:,:,:)     = ZERO*density_hat_real(1,1,1)
        vel2_hat_imag_gq(:,:,:)     = ZERO*density_hat_real(1,1,1)
        vel3_hat_real_gq(:,:,:)     = ZERO*density_hat_real(1,1,1)
        vel3_hat_imag_gq(:,:,:)     = ZERO*density_hat_real(1,1,1)
        pressure_hat_real_gq(:,:,:) = ZERO*density_hat_real(1,1,1)
        pressure_hat_imag_gq(:,:,:) = ZERO*density_hat_real(1,1,1)
        do igq = 1,size(coords)
            do itheta = 1,size(density_hat_real,2)
                do itime = 1,size(density_hat_real,3)
                    density_hat_real_gq( igq,itheta,itime) = interpolate_linear_ad(self%r,density_hat_real(:,itheta,itime), coords(igq)%c1_)
                    density_hat_imag_gq( igq,itheta,itime) = interpolate_linear_ad(self%r,density_hat_imag(:,itheta,itime), coords(igq)%c1_)
                    vel1_hat_real_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel1_hat_real(:,itheta,itime),    coords(igq)%c1_)
                    vel1_hat_imag_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel1_hat_imag(:,itheta,itime),    coords(igq)%c1_)
                    vel2_hat_real_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel2_hat_real(:,itheta,itime),    coords(igq)%c1_)
                    vel2_hat_imag_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel2_hat_imag(:,itheta,itime),    coords(igq)%c1_)
                    vel3_hat_real_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel3_hat_real(:,itheta,itime),    coords(igq)%c1_)
                    vel3_hat_imag_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel3_hat_imag(:,itheta,itime),    coords(igq)%c1_)
                    pressure_hat_real_gq(igq,itheta,itime) = interpolate_linear_ad(self%r,pressure_hat_real(:,itheta,itime),coords(igq)%c1_)
                    pressure_hat_imag_gq(igq,itheta,itime) = interpolate_linear_ad(self%r,pressure_hat_imag(:,itheta,itime),coords(igq)%c1_)
                end do !itime
            end do !itheta
        end do !igq

    end subroutine interpolate_raux_to_rgq
    !********************************************************************************







    !>  DANIEL'S FORMULATION
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine primitive_to_eigenmodes(self,worker,bc_comm,             &
                                       density_bar_r,                   &
                                       vel1_bar_r,                      &
                                       vel2_bar_r,                      &
                                       vel3_bar_r,                      &
                                       pressure_bar_r,                  &
                                       c_bar_r,                         &
                                       density_real,  density_imag,     &
                                       vel1_real,     vel1_imag,        &
                                       vel2_real,     vel2_imag,        &
                                       vel3_real,     vel3_imag,        &
                                       pressure_real, pressure_imag,    &
                                       a1_real,       a1_imag,          &
                                       a2_real,       a2_imag,          &
                                       a3_real,       a3_imag,          &
                                       a4_real,       a4_imag,          &
                                       a5_real,       a5_imag)
        class(nonlocal_nrbc_lindblad_base_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(mpi_comm),             intent(in)      :: bc_comm
        type(AD_D),                 intent(in)      :: density_bar_r(:)
        type(AD_D),                 intent(in)      :: vel1_bar_r(:)
        type(AD_D),                 intent(in)      :: vel2_bar_r(:)
        type(AD_D),                 intent(in)      :: vel3_bar_r(:)
        type(AD_D),                 intent(in)      :: pressure_bar_r(:)
        type(AD_D),                 intent(in)      :: c_bar_r(:)
        type(AD_D),                 intent(inout)   :: density_real(:,:,:)
        type(AD_D),                 intent(inout)   :: density_imag(:,:,:)
        type(AD_D),                 intent(inout)   :: vel1_real(:,:,:)
        type(AD_D),                 intent(inout)   :: vel1_imag(:,:,:)
        type(AD_D),                 intent(inout)   :: vel2_real(:,:,:)
        type(AD_D),                 intent(inout)   :: vel2_imag(:,:,:)
        type(AD_D),                 intent(inout)   :: vel3_real(:,:,:)
        type(AD_D),                 intent(inout)   :: vel3_imag(:,:,:)
        type(AD_D),                 intent(inout)   :: pressure_real(:,:,:)
        type(AD_D),                 intent(inout)   :: pressure_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a1_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a1_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a2_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a2_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a3_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a3_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a4_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a4_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a5_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: a5_imag(:,:,:)

        type(AD_D)  :: k1, k2, k3, k4_real, k4_imag, k5_real, k5_imag, &
                       density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, c_bar, &
                       c_real, c_imag, denom, Tinv_real(5,5), Tinv_imag(5,5), beta

        real(rk),       allocatable, dimension(:)   :: unorm3
        real(rk)                                    :: theta_offset, omega, kz, lm
        integer(ik)                                 :: iradius, igq, ierr, itheta, ntheta, itime, ntime, nr
        logical                                     :: space_time_average

        nr     = size(density_real,1)
        ntheta = size(density_real,2)
        ntime  = size(density_real,3)

        ! Allocate storage
        a1_real = ZERO*density_real
        a2_real = ZERO*density_real
        a3_real = ZERO*density_real
        a4_real = ZERO*density_real
        a5_real = ZERO*density_real
        a1_imag = ZERO*density_real
        a2_imag = ZERO*density_real
        a3_imag = ZERO*density_real
        a4_imag = ZERO*density_real
        a5_imag = ZERO*density_real


        ! Project
        call write_line('WARNING: CHECK DEFINITION OF lm PITCH.', io_proc=IRANK, silence=(verbosity<5))
        do iradius = 1,nr
            ! Get radius-local average
            density_bar  = density_bar_r(iradius)
            vel1_bar     = vel1_bar_r(iradius)
            vel2_bar     = vel2_bar_r(iradius)
            vel3_bar     = vel3_bar_r(iradius)
            pressure_bar = pressure_bar_r(iradius)
            c_bar        = c_bar_r(iradius)

            ! starting with 2 here because the first mode is treated with 1D characteristics
            do itheta = 1,ntheta
                do itime = 1,ntime
                    
                    ! Space-time average handled at the bottom
                    if (itime == 1 .and. itheta == 1) then


                    else


                        ! Get temporal/spatial frequencies
                        omega = get_omega(worker,itime)
                        lm    = get_lm(itheta,ntheta)

                        ! Compute wavenumbers
                        call self%compute_eigenvalues(worker,lm,omega,vel1_bar,vel2_bar,vel3_bar,c_bar, &
                                                      kz, k1, k2, k3, k4_real, k4_imag, k5_real, k5_imag)

                        ! Zero Tinv
                        Tinv_real = ZERO*density_real(1,1,1)
                        Tinv_imag = ZERO*density_real(1,1,1)


                        ! Assemble (1,:)
                        Tinv_real(1,1) = ONE/density_bar
                        Tinv_real(1,5) = -ONE/(density_bar*c_bar*c_bar)

                        ! Assemble (2,:)
                        Tinv_real(2,3) = ONE/c_bar

                        ! Assemble (3,:)
                        Tinv_real(3,2) = (-kz/(c_bar*(k1*k1+kz*kz)))
                        Tinv_real(3,4) = ( k1/(c_bar*(k1*k1+kz*kz)))
                        Tinv_real(3,5) = (-kz/(density_bar*c_bar*vel3_bar*(k1*k1+kz*kz)))


                        ! Assemble (4,:)
                        ! Contribution from vel3:
                        denom  = TWO*c_bar*c_bar*((k4_real*k1 + kz*kz)**TWO + k4_imag*k4_imag*k1*k1)
                        c_real = -vel3_bar*( (k4_real*k1 - k1*k1)*(k4_real*k1 + kz*kz) + (k4_imag*k4_imag*k1*k1) )/denom
                        c_imag = -vel3_bar*( (k4_real*k1 + kz*kz)*k4_imag*k1 - (k4_real*k1 - k1*k1)*k4_imag*k1 )/denom
                        Tinv_real(4,2) = c_real
                        Tinv_imag(4,2) = c_imag

                        ! Contribution from vel2: 
                        c_real = -vel3_bar*( (k4_real*kz - k1*kz)*(k4_real*k1 + kz*kz) + (k4_imag*k4_imag*kz*k1) )/denom
                        c_imag = -vel3_bar*( (k4_real*k1 + kz*kz)*k4_imag*kz - (k4_real*kz - k1*kz)*k4_imag*k1 )/denom
                        Tinv_real(4,4) = c_real
                        Tinv_imag(4,4) = c_imag

                        ! Contribution from pressure:
                        Tinv_real(4,5) = ONE/(TWO*density_bar*c_bar*c_bar)


                        ! Assemble (5,:)
                        ! Contribution from vel3:
                        denom  = TWO*c_bar*c_bar*((k5_real*k1 + kz*kz)**TWO + k5_imag*k5_imag*k1*k1)
                        c_real = -vel3_bar*( (k5_real*k1 - k1*k1)*(k5_real*k1 + kz*kz) + (k5_imag*k5_imag*k1*k1) )/denom
                        c_imag = -vel3_bar*( (k5_real*k1 + kz*kz)*k5_imag*k1 - (k5_real*k1 - k1*k1)*k5_imag*k1 )/denom
                        Tinv_real(5,2) = c_real
                        Tinv_imag(5,2) = c_imag

                        ! Contribution from vel2: 
                        c_real = -vel3_bar*( (k5_real*kz - k1*kz)*(k5_real*k1 + kz*kz) + (k5_imag*k5_imag*kz*k1) )/denom
                        c_imag = -vel3_bar*( (k5_real*k1 + kz*kz)*k5_imag*kz - (k5_real*kz - k1*kz)*k5_imag*k1 )/denom
                        Tinv_real(5,4) = c_real
                        Tinv_imag(5,4) = c_imag

                        ! Contribution from pressure:
                        Tinv_real(5,5) = ONE/(TWO*density_bar*c_bar*c_bar)

                        ! Project
                        a1_real(iradius,itheta,itime) = Tinv_real(1,1)*density_real(iradius,itheta,itime)  - Tinv_imag(1,1)*density_imag(iradius,itheta,itime) + &
                                                        Tinv_real(1,2)*vel3_real(iradius,itheta,itime)     - Tinv_imag(1,2)*vel3_imag(iradius,itheta,itime)    + &
                                                        Tinv_real(1,3)*vel1_real(iradius,itheta,itime)     - Tinv_imag(1,3)*vel1_imag(iradius,itheta,itime)    + &
                                                        Tinv_real(1,4)*vel2_real(iradius,itheta,itime)     - Tinv_imag(1,4)*vel2_imag(iradius,itheta,itime)    + &
                                                        Tinv_real(1,5)*pressure_real(iradius,itheta,itime) - Tinv_imag(1,5)*pressure_imag(iradius,itheta,itime)

                        a1_imag(iradius,itheta,itime) = Tinv_real(1,1)*density_imag(iradius,itheta,itime)  + Tinv_imag(1,1)*density_real(iradius,itheta,itime) + &
                                                        Tinv_real(1,2)*vel3_imag(iradius,itheta,itime)     + Tinv_imag(1,2)*vel3_real(iradius,itheta,itime)    + &
                                                        Tinv_real(1,3)*vel1_imag(iradius,itheta,itime)     + Tinv_imag(1,3)*vel1_real(iradius,itheta,itime)    + &
                                                        Tinv_real(1,4)*vel2_imag(iradius,itheta,itime)     + Tinv_imag(1,4)*vel2_real(iradius,itheta,itime)    + &
                                                        Tinv_real(1,5)*pressure_imag(iradius,itheta,itime) + Tinv_imag(1,5)*pressure_real(iradius,itheta,itime)


                        a2_real(iradius,itheta,itime) = Tinv_real(2,1)*density_real(iradius,itheta,itime)  - Tinv_imag(2,1)*density_imag(iradius,itheta,itime) + &
                                                        Tinv_real(2,2)*vel3_real(iradius,itheta,itime)     - Tinv_imag(2,2)*vel3_imag(iradius,itheta,itime)    + &
                                                        Tinv_real(2,3)*vel1_real(iradius,itheta,itime)     - Tinv_imag(2,3)*vel1_imag(iradius,itheta,itime)    + &
                                                        Tinv_real(2,4)*vel2_real(iradius,itheta,itime)     - Tinv_imag(2,4)*vel2_imag(iradius,itheta,itime)    + &
                                                        Tinv_real(2,5)*pressure_real(iradius,itheta,itime) - Tinv_imag(2,5)*pressure_imag(iradius,itheta,itime)

                        a2_imag(iradius,itheta,itime) = Tinv_real(2,1)*density_imag(iradius,itheta,itime)  + Tinv_imag(2,1)*density_real(iradius,itheta,itime) + &
                                                        Tinv_real(2,2)*vel3_imag(iradius,itheta,itime)     + Tinv_imag(2,2)*vel3_real(iradius,itheta,itime)    + &
                                                        Tinv_real(2,3)*vel1_imag(iradius,itheta,itime)     + Tinv_imag(2,3)*vel1_real(iradius,itheta,itime)    + &
                                                        Tinv_real(2,4)*vel2_imag(iradius,itheta,itime)     + Tinv_imag(2,4)*vel2_real(iradius,itheta,itime)    + &
                                                        Tinv_real(2,5)*pressure_imag(iradius,itheta,itime) + Tinv_imag(2,5)*pressure_real(iradius,itheta,itime)


                        a3_real(iradius,itheta,itime) = Tinv_real(3,1)*density_real(iradius,itheta,itime)  - Tinv_imag(3,1)*density_imag(iradius,itheta,itime) + &
                                                        Tinv_real(3,2)*vel3_real(iradius,itheta,itime)     - Tinv_imag(3,2)*vel3_imag(iradius,itheta,itime)    + &
                                                        Tinv_real(3,3)*vel1_real(iradius,itheta,itime)     - Tinv_imag(3,3)*vel1_imag(iradius,itheta,itime)    + &
                                                        Tinv_real(3,4)*vel2_real(iradius,itheta,itime)     - Tinv_imag(3,4)*vel2_imag(iradius,itheta,itime)    + &
                                                        Tinv_real(3,5)*pressure_real(iradius,itheta,itime) - Tinv_imag(3,5)*pressure_imag(iradius,itheta,itime)

                        a3_imag(iradius,itheta,itime) = Tinv_real(3,1)*density_imag(iradius,itheta,itime)  + Tinv_imag(3,1)*density_real(iradius,itheta,itime) + &
                                                        Tinv_real(3,2)*vel3_imag(iradius,itheta,itime)     + Tinv_imag(3,2)*vel3_real(iradius,itheta,itime)    + &
                                                        Tinv_real(3,3)*vel1_imag(iradius,itheta,itime)     + Tinv_imag(3,3)*vel1_real(iradius,itheta,itime)    + &
                                                        Tinv_real(3,4)*vel2_imag(iradius,itheta,itime)     + Tinv_imag(3,4)*vel2_real(iradius,itheta,itime)    + &
                                                        Tinv_real(3,5)*pressure_imag(iradius,itheta,itime) + Tinv_imag(3,5)*pressure_real(iradius,itheta,itime)


                        a4_real(iradius,itheta,itime) = Tinv_real(4,1)*density_real(iradius,itheta,itime)  - Tinv_imag(4,1)*density_imag(iradius,itheta,itime) + &
                                                        Tinv_real(4,2)*vel3_real(iradius,itheta,itime)     - Tinv_imag(4,2)*vel3_imag(iradius,itheta,itime)    + &
                                                        Tinv_real(4,3)*vel1_real(iradius,itheta,itime)     - Tinv_imag(4,3)*vel1_imag(iradius,itheta,itime)    + &
                                                        Tinv_real(4,4)*vel2_real(iradius,itheta,itime)     - Tinv_imag(4,4)*vel2_imag(iradius,itheta,itime)    + &
                                                        Tinv_real(4,5)*pressure_real(iradius,itheta,itime) - Tinv_imag(4,5)*pressure_imag(iradius,itheta,itime)

                        a4_imag(iradius,itheta,itime) = Tinv_real(4,1)*density_imag(iradius,itheta,itime)  + Tinv_imag(4,1)*density_real(iradius,itheta,itime) + &
                                                        Tinv_real(4,2)*vel3_imag(iradius,itheta,itime)     + Tinv_imag(4,2)*vel3_real(iradius,itheta,itime)    + &
                                                        Tinv_real(4,3)*vel1_imag(iradius,itheta,itime)     + Tinv_imag(4,3)*vel1_real(iradius,itheta,itime)    + &
                                                        Tinv_real(4,4)*vel2_imag(iradius,itheta,itime)     + Tinv_imag(4,4)*vel2_real(iradius,itheta,itime)    + &
                                                        Tinv_real(4,5)*pressure_imag(iradius,itheta,itime) + Tinv_imag(4,5)*pressure_real(iradius,itheta,itime)


                        a5_real(iradius,itheta,itime) = Tinv_real(5,1)*density_real(iradius,itheta,itime)  - Tinv_imag(5,1)*density_imag(iradius,itheta,itime) + &
                                                        Tinv_real(5,2)*vel3_real(iradius,itheta,itime)     - Tinv_imag(5,2)*vel3_imag(iradius,itheta,itime)    + &
                                                        Tinv_real(5,3)*vel1_real(iradius,itheta,itime)     - Tinv_imag(5,3)*vel1_imag(iradius,itheta,itime)    + &
                                                        Tinv_real(5,4)*vel2_real(iradius,itheta,itime)     - Tinv_imag(5,4)*vel2_imag(iradius,itheta,itime)    + &
                                                        Tinv_real(5,5)*pressure_real(iradius,itheta,itime) - Tinv_imag(5,5)*pressure_imag(iradius,itheta,itime)

                        a5_imag(iradius,itheta,itime) = Tinv_real(5,1)*density_imag(iradius,itheta,itime)  + Tinv_imag(5,1)*density_real(iradius,itheta,itime) + &
                                                        Tinv_real(5,2)*vel3_imag(iradius,itheta,itime)     + Tinv_imag(5,2)*vel3_real(iradius,itheta,itime)    + &
                                                        Tinv_real(5,3)*vel1_imag(iradius,itheta,itime)     + Tinv_imag(5,3)*vel1_real(iradius,itheta,itime)    + &
                                                        Tinv_real(5,4)*vel2_imag(iradius,itheta,itime)     + Tinv_imag(5,4)*vel2_real(iradius,itheta,itime)    + &
                                                        Tinv_real(5,5)*pressure_imag(iradius,itheta,itime) + Tinv_imag(5,5)*pressure_real(iradius,itheta,itime)

                    end if
                end do !itime
            end do !itheta
        end do !iradius



    end subroutine primitive_to_eigenmodes
    !********************************************************************************





    !>  DANIEL'S FORMULATION
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine eigenmodes_to_primitive(self,worker,bc_comm,             &
                                       density_bar_r, vel1_bar_r, vel2_bar_r, vel3_bar_r, pressure_bar_r, c_bar_r, &
                                       a1_real,       a1_imag,          &
                                       a2_real,       a2_imag,          &
                                       a3_real,       a3_imag,          &
                                       a4_real,       a4_imag,          &
                                       a5_real,       a5_imag,          &
                                       density_real,  density_imag,     &
                                       vel1_real,     vel1_imag,        &
                                       vel2_real,     vel2_imag,        &
                                       vel3_real,     vel3_imag,        &
                                       pressure_real, pressure_imag)
        class(nonlocal_nrbc_lindblad_base_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(mpi_comm),             intent(in)      :: bc_comm
        type(AD_D),                 intent(in)      :: density_bar_r(:)
        type(AD_D),                 intent(in)      :: vel1_bar_r(:)
        type(AD_D),                 intent(in)      :: vel2_bar_r(:)
        type(AD_D),                 intent(in)      :: vel3_bar_r(:)
        type(AD_D),                 intent(in)      :: pressure_bar_r(:)
        type(AD_D),                 intent(in)      :: c_bar_r(:)
        type(AD_D),                 intent(in)      :: a1_real(:,:,:)
        type(AD_D),                 intent(in)      :: a1_imag(:,:,:)
        type(AD_D),                 intent(in)      :: a2_real(:,:,:)
        type(AD_D),                 intent(in)      :: a2_imag(:,:,:)
        type(AD_D),                 intent(in)      :: a3_real(:,:,:)
        type(AD_D),                 intent(in)      :: a3_imag(:,:,:)
        type(AD_D),                 intent(in)      :: a4_real(:,:,:)
        type(AD_D),                 intent(in)      :: a4_imag(:,:,:)
        type(AD_D),                 intent(in)      :: a5_real(:,:,:)
        type(AD_D),                 intent(in)      :: a5_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: density_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: density_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_imag(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_real(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_imag(:,:,:)


        type(AD_D)  :: k1, k2, k3, k4_real, k4_imag, k5_real, k5_imag, &
                       density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, &
                       c_bar,  denom, c_real, c_imag, T_real(5,5), T_imag(5,5), beta
        
        real(rk),       allocatable, dimension(:)   :: unorm3
        real(rk)                                    :: theta_offset, omega, kz, lm
        integer(ik)                                 :: iradius, igq, ierr, itheta, ntheta, itime, ntime, nr
        logical                                     :: space_time_average


        density_real  = ZERO*a1_real
        vel1_real     = ZERO*a1_real
        vel2_real     = ZERO*a1_real
        vel3_real     = ZERO*a1_real
        pressure_real = ZERO*a1_real
        density_imag  = ZERO*a1_real
        vel1_imag     = ZERO*a1_real
        vel2_imag     = ZERO*a1_real
        vel3_imag     = ZERO*a1_real
        pressure_imag = ZERO*a1_real


        nr     = size(density_real,1)
        ntheta = size(density_real,2)
        ntime  = size(density_real,3)

        ! Project
        call write_line('WARNING: CHECK DEFINITION OF lm PITCH.', io_proc=IRANK, silence=(verbosity<5))
        do iradius = 1,nr
            ! Get radius-local average
            density_bar  = density_bar_r(iradius)
            vel1_bar     = vel1_bar_r(iradius)
            vel2_bar     = vel2_bar_r(iradius)
            vel3_bar     = vel3_bar_r(iradius)
            pressure_bar = pressure_bar_r(iradius)
            c_bar        = c_bar_r(iradius)

            ! starting with 2 here because the first mode is treated with 1D characteristics
            do itheta = 1,ntheta
                do itime = 1,ntime
                    
                    ! Space-time average handled at the bottom
                    if (itime == 1 .and. itheta == 1) then


                    else

                        ! Get temporal/spatial frequencies
                        omega = get_omega(worker,itime)
                        lm    = get_lm(itheta,ntheta)

                        ! Compute wavenumbers
                        call self%compute_eigenvalues(worker,lm,omega,vel1_bar,vel2_bar,vel3_bar,c_bar, &
                                                      kz, k1, k2, k3, k4_real, k4_imag, k5_real, k5_imag)

                        ! First zero fields
                        density_real(iradius,itheta,itime)  = ZERO
                        density_imag(iradius,itheta,itime)  = ZERO
                        vel1_real(iradius,itheta,itime)     = ZERO
                        vel1_imag(iradius,itheta,itime)     = ZERO
                        vel2_real(iradius,itheta,itime)     = ZERO
                        vel2_imag(iradius,itheta,itime)     = ZERO
                        vel3_real(iradius,itheta,itime)     = ZERO
                        vel3_imag(iradius,itheta,itime)     = ZERO
                        pressure_real(iradius,itheta,itime) = ZERO
                        pressure_imag(iradius,itheta,itime) = ZERO
                        T_real = ZERO*density_real(1,1,1)
                        T_imag = ZERO*density_real(1,1,1)


                        ! Assemble (1,:)
                        T_real(1,1) = density_bar
                        T_real(1,4) = density_bar
                        T_real(1,5) = density_bar

                        
                        ! Assemble (2,:)
                        T_real(2,3) = -c_bar*kz

                        ! Contribution from a4
                        denom = vel3_bar*( (k4_real-k1)**TWO + k4_imag*k4_imag)
                        c_real = (-c_bar*c_bar*(k4_real*(k4_real-k1) + k4_imag*k4_imag)/denom)
                        c_imag = (c_bar*c_bar*(k4_imag*k1)/denom)
                        T_real(2,4) = c_real
                        T_imag(2,4) = c_imag

                        ! Contribution from a5
                        denom = vel3_bar*( (k5_real-k1)**TWO + k5_imag*k5_imag)
                        c_real = (-c_bar*c_bar*(k5_real*(k5_real-k1) + k5_imag*k5_imag)/denom)
                        c_imag = (c_bar*c_bar*(k5_imag*k1)/denom)
                        T_real(2,5) = c_real
                        T_imag(2,5) = c_imag


                        ! Assemble (3,:)
                        T_real(3,2) = c_bar


                        ! Assemble (4,:)
                        ! Contribution from a3
                        T_real(4,3) = c_bar*k1

                        ! Contribution from a4
                        denom = vel3_bar*((k4_real-k1)**TWO + k4_imag*k4_imag)
                        c_real = -c_bar*c_bar*kz*(k4_real-k1)/denom
                        c_imag =  c_bar*c_bar*kz*k4_imag/denom
                        T_real(4,4) = c_real
                        T_imag(4,4) = c_imag

                        ! Contribution from a5
                        denom = vel3_bar*((k5_real-k1)**TWO + k5_imag*k5_imag)
                        c_real = -c_bar*c_bar*kz*(k5_real-k1)/denom
                        c_imag =  c_bar*c_bar*kz*k5_imag/denom
                        T_real(4,5) = c_real
                        T_imag(4,5) = c_imag


                        ! Assemble (5,:)
                        T_real(5,4) = density_bar*c_bar*c_bar
                        T_real(5,5) = density_bar*c_bar*c_bar


                        ! Density
                        density_real(iradius,itheta,itime) = T_real(1,1)*a1_real(iradius,itheta,itime) - T_imag(1,1)*a1_imag(iradius,itheta,itime) + &
                                                             T_real(1,2)*a2_real(iradius,itheta,itime) - T_imag(1,2)*a2_imag(iradius,itheta,itime) + &
                                                             T_real(1,3)*a3_real(iradius,itheta,itime) - T_imag(1,3)*a3_imag(iradius,itheta,itime) + &
                                                             T_real(1,4)*a4_real(iradius,itheta,itime) - T_imag(1,4)*a4_imag(iradius,itheta,itime) + &
                                                             T_real(1,5)*a5_real(iradius,itheta,itime) - T_imag(1,5)*a5_imag(iradius,itheta,itime)

                        density_imag(iradius,itheta,itime) = T_real(1,1)*a1_imag(iradius,itheta,itime) + T_imag(1,1)*a1_real(iradius,itheta,itime) + &
                                                             T_real(1,2)*a2_imag(iradius,itheta,itime) + T_imag(1,2)*a2_real(iradius,itheta,itime) + &
                                                             T_real(1,3)*a3_imag(iradius,itheta,itime) + T_imag(1,3)*a3_real(iradius,itheta,itime) + &
                                                             T_real(1,4)*a4_imag(iradius,itheta,itime) + T_imag(1,4)*a4_real(iradius,itheta,itime) + &
                                                             T_real(1,5)*a5_imag(iradius,itheta,itime) + T_imag(1,5)*a5_real(iradius,itheta,itime)

                        ! Vel3
                        vel3_real(iradius,itheta,itime) = T_real(2,1)*a1_real(iradius,itheta,itime) - T_imag(2,1)*a1_imag(iradius,itheta,itime) + &
                                                          T_real(2,2)*a2_real(iradius,itheta,itime) - T_imag(2,2)*a2_imag(iradius,itheta,itime) + &
                                                          T_real(2,3)*a3_real(iradius,itheta,itime) - T_imag(2,3)*a3_imag(iradius,itheta,itime) + &
                                                          T_real(2,4)*a4_real(iradius,itheta,itime) - T_imag(2,4)*a4_imag(iradius,itheta,itime) + &
                                                          T_real(2,5)*a5_real(iradius,itheta,itime) - T_imag(2,5)*a5_imag(iradius,itheta,itime)

                        vel3_imag(iradius,itheta,itime) = T_real(2,1)*a1_imag(iradius,itheta,itime) + T_imag(2,1)*a1_real(iradius,itheta,itime) + &
                                                          T_real(2,2)*a2_imag(iradius,itheta,itime) + T_imag(2,2)*a2_real(iradius,itheta,itime) + &
                                                          T_real(2,3)*a3_imag(iradius,itheta,itime) + T_imag(2,3)*a3_real(iradius,itheta,itime) + &
                                                          T_real(2,4)*a4_imag(iradius,itheta,itime) + T_imag(2,4)*a4_real(iradius,itheta,itime) + &
                                                          T_real(2,5)*a5_imag(iradius,itheta,itime) + T_imag(2,5)*a5_real(iradius,itheta,itime)

                        ! Vel1
                        vel1_real(iradius,itheta,itime) = T_real(3,1)*a1_real(iradius,itheta,itime) - T_imag(3,1)*a1_imag(iradius,itheta,itime) + &
                                                          T_real(3,2)*a2_real(iradius,itheta,itime) - T_imag(3,2)*a2_imag(iradius,itheta,itime) + &
                                                          T_real(3,3)*a3_real(iradius,itheta,itime) - T_imag(3,3)*a3_imag(iradius,itheta,itime) + &
                                                          T_real(3,4)*a4_real(iradius,itheta,itime) - T_imag(3,4)*a4_imag(iradius,itheta,itime) + &
                                                          T_real(3,5)*a5_real(iradius,itheta,itime) - T_imag(3,5)*a5_imag(iradius,itheta,itime)

                        vel1_imag(iradius,itheta,itime) = T_real(3,1)*a1_imag(iradius,itheta,itime) + T_imag(3,1)*a1_real(iradius,itheta,itime) + &
                                                          T_real(3,2)*a2_imag(iradius,itheta,itime) + T_imag(3,2)*a2_real(iradius,itheta,itime) + &
                                                          T_real(3,3)*a3_imag(iradius,itheta,itime) + T_imag(3,3)*a3_real(iradius,itheta,itime) + &
                                                          T_real(3,4)*a4_imag(iradius,itheta,itime) + T_imag(3,4)*a4_real(iradius,itheta,itime) + &
                                                          T_real(3,5)*a5_imag(iradius,itheta,itime) + T_imag(3,5)*a5_real(iradius,itheta,itime)

                        ! Vel2
                        vel2_real(iradius,itheta,itime) = T_real(4,1)*a1_real(iradius,itheta,itime) - T_imag(4,1)*a1_imag(iradius,itheta,itime) + &
                                                          T_real(4,2)*a2_real(iradius,itheta,itime) - T_imag(4,2)*a2_imag(iradius,itheta,itime) + &
                                                          T_real(4,3)*a3_real(iradius,itheta,itime) - T_imag(4,3)*a3_imag(iradius,itheta,itime) + &
                                                          T_real(4,4)*a4_real(iradius,itheta,itime) - T_imag(4,4)*a4_imag(iradius,itheta,itime) + &
                                                          T_real(4,5)*a5_real(iradius,itheta,itime) - T_imag(4,5)*a5_imag(iradius,itheta,itime)

                        vel2_imag(iradius,itheta,itime) = T_real(4,1)*a1_imag(iradius,itheta,itime) + T_imag(4,1)*a1_real(iradius,itheta,itime) + &
                                                          T_real(4,2)*a2_imag(iradius,itheta,itime) + T_imag(4,2)*a2_real(iradius,itheta,itime) + &
                                                          T_real(4,3)*a3_imag(iradius,itheta,itime) + T_imag(4,3)*a3_real(iradius,itheta,itime) + &
                                                          T_real(4,4)*a4_imag(iradius,itheta,itime) + T_imag(4,4)*a4_real(iradius,itheta,itime) + &
                                                          T_real(4,5)*a5_imag(iradius,itheta,itime) + T_imag(4,5)*a5_real(iradius,itheta,itime)

                        ! Pressure
                        pressure_real(iradius,itheta,itime) = T_real(5,1)*a1_real(iradius,itheta,itime) - T_imag(5,1)*a1_imag(iradius,itheta,itime) + &
                                                              T_real(5,2)*a2_real(iradius,itheta,itime) - T_imag(5,2)*a2_imag(iradius,itheta,itime) + &
                                                              T_real(5,3)*a3_real(iradius,itheta,itime) - T_imag(5,3)*a3_imag(iradius,itheta,itime) + &
                                                              T_real(5,4)*a4_real(iradius,itheta,itime) - T_imag(5,4)*a4_imag(iradius,itheta,itime) + &
                                                              T_real(5,5)*a5_real(iradius,itheta,itime) - T_imag(5,5)*a5_imag(iradius,itheta,itime)

                        pressure_imag(iradius,itheta,itime) = T_real(5,1)*a1_imag(iradius,itheta,itime) + T_imag(5,1)*a1_real(iradius,itheta,itime) + &
                                                              T_real(5,2)*a2_imag(iradius,itheta,itime) + T_imag(5,2)*a2_real(iradius,itheta,itime) + &
                                                              T_real(5,3)*a3_imag(iradius,itheta,itime) + T_imag(5,3)*a3_real(iradius,itheta,itime) + &
                                                              T_real(5,4)*a4_imag(iradius,itheta,itime) + T_imag(5,4)*a4_real(iradius,itheta,itime) + &
                                                              T_real(5,5)*a5_imag(iradius,itheta,itime) + T_imag(5,5)*a5_real(iradius,itheta,itime)

                    end if

                end do !itime
            end do !itheta
        end do !iradius


    end subroutine eigenmodes_to_primitive
    !********************************************************************************
    



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/8/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_eigenvalues(self, worker, lm, omega, vel1_bar, vel2_bar, vel3_bar, c_bar, &
                                   kz, k1, k2, k3, k4_real, k4_imag, k5_real, k5_imag)
        class(nonlocal_nrbc_lindblad_base_t), intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        real(rk),               intent(in)      :: lm
        real(rk),               intent(in)      :: omega
        type(AD_D),             intent(in)      :: vel1_bar
        type(AD_D),             intent(in)      :: vel2_bar
        type(AD_D),             intent(in)      :: vel3_bar
        type(AD_D),             intent(in)      :: c_bar
        real(rk),               intent(inout)   :: kz
        type(AD_D),             intent(inout)   :: k1
        type(AD_D),             intent(inout)   :: k2
        type(AD_D),             intent(inout)   :: k3
        type(AD_D),             intent(inout)   :: k4_real
        type(AD_D),             intent(inout)   :: k4_imag
        type(AD_D),             intent(inout)   :: k5_real
        type(AD_D),             intent(inout)   :: k5_imag

        real(rk),   allocatable, dimension(:)   :: pitch, unorm3
        type(AD_D)  :: pyra, ktmp_real, ktmp_imag, vg4, vg5

        complex(rk) :: omega_c, pyra_c, k4_c, k5_c
        real(rk)    :: vel2_bar_r, vel3_bar_r, c_bar_r

        pitch  = self%bcproperties%compute('Pitch A',worker%time(),worker%coords())

        kz = TWO*PI*lm/pitch(1)
        pyra = (omega - kz*vel2_bar)**TWO - kz*kz*(c_bar**TWO - vel3_bar**TWO)

        ! Compute k1,k2,k3
        k1 = (omega - kz*vel2_bar)/vel3_bar
        k2 = (omega - kz*vel2_bar)/vel3_bar
        k3 = (omega - kz*vel2_bar)/vel3_bar

        if (pyra >= ZERO) then
            k4_real = (-vel3_bar*(omega - kz*vel2_bar) + c_bar*sqrt(pyra))/(c_bar**TWO - vel3_bar**TWO)
            k5_real = (-vel3_bar*(omega - kz*vel2_bar) - c_bar*sqrt(pyra))/(c_bar**TWO - vel3_bar**TWO)
            k4_imag = ZERO*k4_real
            k5_imag = ZERO*k4_real

!            ! test direction of propagation by perturbing omega
!            vel3_bar_r = vel3_bar%x_ad_
!            vel2_bar_r = vel2_bar%x_ad_
!            c_bar_r = c_bar%x_ad_
!            omega_c = cmplx(omega,-0.00001_rk)
!            pyra_c = (omega_c - kz*vel2_bar_r)**TWO - kz*kz*(c_bar_r**TWO - vel3_bar_r**TWO)
!            k4_c = (-vel3_bar_r*(omega_c - kz*vel2_bar_r) + c_bar_r*sqrt(pyra_c))/(c_bar_r**TWO - vel3_bar_r**TWO)
!            k5_c = (-vel3_bar_r*(omega_c - kz*vel2_bar_r) - c_bar_r*sqrt(pyra_c))/(c_bar_r**TWO - vel3_bar_r**TWO)
!
!            if (imagpart(k4_c) < ZERO) then
!                k4_imag = -0.0000001_rk
!            else
!                k4_imag = 0.0000001_rk
!            end if
!
!            if (imagpart(k5_c) < ZERO) then
!                k5_imag = -0.0000001_rk
!            else
!                k5_imag = 0.0000001_rk
!            end if


            vg4 = -(c_bar*c_bar - vel3_bar*vel3_bar)/(vel3_bar - (c_bar*(omega-kz*vel2_bar)/sqrt(pyra)))
            vg5 = -(c_bar*c_bar - vel3_bar*vel3_bar)/(vel3_bar + (c_bar*(omega-kz*vel2_bar)/sqrt(pyra)))


        else
            k4_real = (-vel3_bar*(omega - kz*vel2_bar))/(c_bar**TWO - vel3_bar**TWO)
            k4_imag =  c_bar*sqrt(-pyra)/(c_bar**TWO - vel3_bar**TWO)
            k5_real = (-vel3_bar*(omega - kz*vel2_bar))/(c_bar**TWO - vel3_bar**TWO)
            k5_imag = -c_bar*sqrt(-pyra)/(c_bar**TWO - vel3_bar**TWO)

!            ! test direction of propagation by perturbing omega
!            vel3_bar_r = vel3_bar%x_ad_
!            vel2_bar_r = vel2_bar%x_ad_
!            c_bar_r = c_bar%x_ad_
!            omega_c = cmplx(omega,-0.00001_rk)
!            pyra_c = (omega_c - kz*vel2_bar_r)**TWO - kz*kz*(c_bar_r**TWO - vel3_bar_r**TWO)
!            k4_c = (-vel3_bar_r*(omega_c - kz*vel2_bar_r) + c_bar_r*sqrt(pyra_c))/(c_bar_r**TWO - vel3_bar_r**TWO)
!            k5_c = (-vel3_bar_r*(omega_c - kz*vel2_bar_r) - c_bar_r*sqrt(pyra_c))/(c_bar_r**TWO - vel3_bar_r**TWO)
!
!            if ( ((imagpart(k4_c) < ZERO) .and. (k4_imag > ZERO)) .or. &
!                 ((imagpart(k4_c) > ZERO) .and. (k4_imag < ZERO)) ) then
!                k4_imag = -k4_imag
!            end if
!
!            if ( ((imagpart(k5_c) < ZERO) .and. (k5_imag > ZERO)) .or. &
!                 ((imagpart(k5_c) > ZERO) .and. (k5_imag < ZERO)) ) then
!                k5_imag = -k5_imag
!            end if

            vg4 = ZERO*k4_real
            vg5 = ZERO*k4_real
            if (k4_imag < ZERO) then
                vg4 = ONE
            else
                vg4 = -ONE
            end if

            if (k5_imag < ZERO) then
                vg5 = ONE
            else
                vg5 = -ONE
            end if

        end if

!        unorm3 = worker%unit_normal(3)
!        if (k5_imag*unorm3(1) < ZERO) then
!            ktmp_real = k5_real
!            ktmp_imag = k5_imag
!
!            k5_real = k4_real
!            k5_imag = k4_imag
!            k4_real = ktmp_real
!            k4_imag = ktmp_imag
!        end if


        unorm3 = worker%unit_normal(3)
        if (unorm3(1) < ZERO) then
            if (vg5 < ZERO) then
                ktmp_real = k5_real
                ktmp_imag = k5_imag

                k5_real = k4_real
                k5_imag = k4_imag
                k4_real = ktmp_real
                k4_imag = ktmp_imag
            end if
        end if

        
        if (unorm3(1) > ZERO) then
            if (vg5 > ZERO) then
                ktmp_real = k5_real
                ktmp_imag = k5_imag

                k5_real = k4_real
                k5_imag = k4_imag
                k4_real = ktmp_real
                k4_imag = ktmp_imag
            end if
        end if

    end subroutine compute_eigenvalues
    !********************************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine primitive_to_characteristics(self,worker,bc_comm,                    &
                                            density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, &
                                            density_hat_real,  density_hat_imag,    &
                                            vel1_hat_real,     vel1_hat_imag,       &
                                            vel2_hat_real,     vel2_hat_imag,       &
                                            vel3_hat_real,     vel3_hat_imag,       &
                                            pressure_hat_real, pressure_hat_imag,   &
                                            c1_hat_real,       c1_hat_imag,         &
                                            c2_hat_real,       c2_hat_imag,         &
                                            c3_hat_real,       c3_hat_imag,         &
                                            c4_hat_real,       c4_hat_imag,         &
                                            c5_hat_real,       c5_hat_imag)
        class(nonlocal_nrbc_lindblad_base_t),         intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        type(mpi_comm),                 intent(in)      :: bc_comm
        type(AD_D),                     intent(in)      :: density_bar(:)
        type(AD_D),                     intent(in)      :: vel1_bar(:)
        type(AD_D),                     intent(in)      :: vel2_bar(:)
        type(AD_D),                     intent(in)      :: vel3_bar(:)
        type(AD_D),                     intent(in)      :: pressure_bar(:)
        type(AD_D),                     intent(in)      :: density_hat_real(:,:,:)
        type(AD_D),                     intent(in)      :: density_hat_imag(:,:,:)
        type(AD_D),                     intent(in)      :: vel1_hat_real(:,:,:)
        type(AD_D),                     intent(in)      :: vel1_hat_imag(:,:,:)
        type(AD_D),                     intent(in)      :: vel2_hat_real(:,:,:)
        type(AD_D),                     intent(in)      :: vel2_hat_imag(:,:,:)
        type(AD_D),                     intent(in)      :: vel3_hat_real(:,:,:)
        type(AD_D),                     intent(in)      :: vel3_hat_imag(:,:,:)
        type(AD_D),                     intent(in)      :: pressure_hat_real(:,:,:)
        type(AD_D),                     intent(in)      :: pressure_hat_imag(:,:,:)
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

        type(AD_D)  :: c_bar

        integer(ik) :: nradius, iradius, itheta, ntheta, ierr, ntime, itime

        ! Define Fourier discretization
        nradius = size(density_hat_real,1)
        ntheta  = size(density_hat_real,2)
        ntime   = size(density_hat_real,3)

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
                    c1_hat_real(iradius,itheta,itime) = -(c_bar*c_bar)*density_hat_real(iradius,itheta,itime)             +  (ONE)*pressure_hat_real(iradius,itheta,itime)
                    c2_hat_real(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel1_hat_real(iradius,itheta,itime)
                    c3_hat_real(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel2_hat_real(iradius,itheta,itime)
                    c4_hat_real(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel3_hat_real(iradius,itheta,itime)  +  (ONE)*pressure_hat_real(iradius,itheta,itime)
                    c5_hat_real(iradius,itheta,itime) = -(density_bar(iradius)*c_bar)*vel3_hat_real(iradius,itheta,itime) +  (ONE)*pressure_hat_real(iradius,itheta,itime)
                                               
                    c1_hat_imag(iradius,itheta,itime) = -(c_bar*c_bar)*density_hat_imag(iradius,itheta,itime)             +  (ONE)*pressure_hat_imag(iradius,itheta,itime)
                    c2_hat_imag(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel1_hat_imag(iradius,itheta,itime)
                    c3_hat_imag(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel2_hat_imag(iradius,itheta,itime)
                    c4_hat_imag(iradius,itheta,itime) = (density_bar(iradius)*c_bar)*vel3_hat_imag(iradius,itheta,itime)  +  (ONE)*pressure_hat_imag(iradius,itheta,itime)
                    c5_hat_imag(iradius,itheta,itime) = -(density_bar(iradius)*c_bar)*vel3_hat_imag(iradius,itheta,itime) +  (ONE)*pressure_hat_imag(iradius,itheta,itime)
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
                                            density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, &
                                            c1_hat_real,       c1_hat_imag,       &
                                            c2_hat_real,       c2_hat_imag,       &
                                            c3_hat_real,       c3_hat_imag,       &
                                            c4_hat_real,       c4_hat_imag,       &
                                            c5_hat_real,       c5_hat_imag,       &
                                            density_hat_real,  density_hat_imag,  &
                                            vel1_hat_real,     vel1_hat_imag,     &
                                            vel2_hat_real,     vel2_hat_imag,     &
                                            vel3_hat_real,     vel3_hat_imag,     &
                                            pressure_hat_real, pressure_hat_imag)
        class(nonlocal_nrbc_lindblad_base_t),         intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        type(mpi_comm),                 intent(in)      :: bc_comm
        type(AD_D),                     intent(in)      :: density_bar(:)
        type(AD_D),                     intent(in)      :: vel1_bar(:)
        type(AD_D),                     intent(in)      :: vel2_bar(:)
        type(AD_D),                     intent(in)      :: vel3_bar(:)
        type(AD_D),                     intent(in)      :: pressure_bar(:)
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
        type(AD_D), allocatable,        intent(inout)   :: density_hat_real(:,:,:)
        type(AD_D), allocatable,        intent(inout)   :: density_hat_imag(:,:,:)
        type(AD_D), allocatable,        intent(inout)   :: vel1_hat_real(:,:,:)
        type(AD_D), allocatable,        intent(inout)   :: vel1_hat_imag(:,:,:)
        type(AD_D), allocatable,        intent(inout)   :: vel2_hat_real(:,:,:)
        type(AD_D), allocatable,        intent(inout)   :: vel2_hat_imag(:,:,:)
        type(AD_D), allocatable,        intent(inout)   :: vel3_hat_real(:,:,:)
        type(AD_D), allocatable,        intent(inout)   :: vel3_hat_imag(:,:,:)
        type(AD_D), allocatable,        intent(inout)   :: pressure_hat_real(:,:,:)
        type(AD_D), allocatable,        intent(inout)   :: pressure_hat_imag(:,:,:)

        type(AD_D)  :: density_bar_r, pressure_bar_r, c_bar_r

        integer(ik) :: iradius, itheta, itime

        density_hat_real  = ZERO*c1_hat_real
        vel1_hat_real     = ZERO*c1_hat_real
        vel2_hat_real     = ZERO*c1_hat_real
        vel3_hat_real     = ZERO*c1_hat_real
        pressure_hat_real = ZERO*c1_hat_real

        density_hat_imag  = ZERO*c1_hat_real
        vel1_hat_imag     = ZERO*c1_hat_real
        vel2_hat_imag     = ZERO*c1_hat_real
        vel3_hat_imag     = ZERO*c1_hat_real
        pressure_hat_imag = ZERO*c1_hat_real


        ! Convert characteristic Fourier modes back to primitive Fourier modes and store
        do iradius = 1,size(c1_hat_real,1)
            ! Get radius-local averages
            density_bar_r  = density_bar(iradius)
            pressure_bar_r = pressure_bar(iradius)
            c_bar_r        = sqrt(gam*pressure_bar_r/density_bar_r)
            do itheta = 1,size(c1_hat_real,2)
                do itime = 1,size(c1_hat_real,3)
                    density_hat_real(iradius,itheta,itime) = -(ONE/(c_bar_r*c_bar_r))*c1_hat_real(iradius,itheta,itime) + (ONE/(TWO*c_bar_r*c_bar_r))*c4_hat_real(iradius,itheta,itime) + (ONE/(TWO*c_bar_r*c_bar_r))*c5_hat_real(iradius,itheta,itime)
                    vel1_hat_real(iradius,itheta,itime) = (ONE/(density_bar_r*c_bar_r))*c2_hat_real(iradius,itheta,itime)
                    vel2_hat_real(iradius,itheta,itime) = (ONE/(density_bar_r*c_bar_r))*c3_hat_real(iradius,itheta,itime)
                    vel3_hat_real(iradius,itheta,itime) = (ONE/(TWO*density_bar_r*c_bar_r))*c4_hat_real(iradius,itheta,itime) - (ONE/(TWO*density_bar_r*c_bar_r))*c5_hat_real(iradius,itheta,itime)
                    pressure_hat_real(iradius,itheta,itime) = HALF*c4_hat_real(iradius,itheta,itime) + HALF*c5_hat_real(iradius,itheta,itime)

                    density_hat_imag(iradius,itheta,itime) = -(ONE/(c_bar_r*c_bar_r))*c1_hat_imag(iradius,itheta,itime) + (ONE/(TWO*c_bar_r*c_bar_r))*c4_hat_imag(iradius,itheta,itime) + (ONE/(TWO*c_bar_r*c_bar_r))*c5_hat_imag(iradius,itheta,itime)
                    vel1_hat_imag(iradius,itheta,itime) = (ONE/(density_bar_r*c_bar_r))*c2_hat_imag(iradius,itheta,itime)
                    vel2_hat_imag(iradius,itheta,itime) = (ONE/(density_bar_r*c_bar_r))*c3_hat_imag(iradius,itheta,itime)
                    vel3_hat_imag(iradius,itheta,itime) = (ONE/(TWO*density_bar_r*c_bar_r))*c4_hat_imag(iradius,itheta,itime) - (ONE/(TWO*density_bar_r*c_bar_r))*c5_hat_imag(iradius,itheta,itime)
                    pressure_hat_imag(iradius,itheta,itime) = HALF*c4_hat_imag(iradius,itheta,itime) + HALF*c5_hat_imag(iradius,itheta,itime)
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
    subroutine compute_boundary_average(self,worker,bc_comm,density_bar,vel1_bar,vel2_bar,vel3_bar,pressure_bar,c_bar, &
                                                            density_avg,vel1_avg,vel2_avg,vel3_avg,pressure_avg,c_avg)
        class(nonlocal_nrbc_lindblad_base_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(in)      :: worker
        type(mpi_comm),                             intent(in)      :: bc_comm
        type(AD_D),                                 intent(in)      :: density_bar(:)
        type(AD_D),                                 intent(in)      :: vel1_bar(:)
        type(AD_D),                                 intent(in)      :: vel2_bar(:)
        type(AD_D),                                 intent(in)      :: vel3_bar(:)
        type(AD_D),                                 intent(in)      :: pressure_bar(:)
        type(AD_D),                                 intent(in)      :: c_bar(:)
        type(AD_D),                                 intent(inout)   :: density_avg
        type(AD_D),                                 intent(inout)   :: vel1_avg
        type(AD_D),                                 intent(inout)   :: vel2_avg
        type(AD_D),                                 intent(inout)   :: vel3_avg
        type(AD_D),                                 intent(inout)   :: pressure_avg
        type(AD_D),                                 intent(inout)   :: c_avg

        real(rk)    :: area, dr
        integer(ik) :: irad


        density_avg  = ZERO*density_bar(1)
        vel1_avg     = ZERO*density_bar(1)
        vel2_avg     = ZERO*density_bar(1)
        vel3_avg     = ZERO*density_bar(1)
        pressure_avg = ZERO*density_bar(1)
        c_avg        = ZERO*density_bar(1)

        if (worker%coordinate_system() == 'Cartesian') then
            dr = self%r(2) - self%r(1)
            area = ZERO
            do irad = 1,size(density_bar)-1
                density_avg  = density_avg  + dr*(density_bar(irad+1)  + density_bar(irad))/TWO
                vel1_avg     = vel1_avg     + dr*(vel1_bar(irad+1)     + vel1_bar(irad))/TWO
                vel2_avg     = vel2_avg     + dr*(vel2_bar(irad+1)     + vel2_bar(irad))/TWO
                vel3_avg     = vel3_avg     + dr*(vel3_bar(irad+1)     + vel3_bar(irad))/TWO
                pressure_avg = pressure_avg + dr*(pressure_bar(irad+1) + pressure_bar(irad))/TWO
                c_avg        = c_avg        + dr*(c_bar(irad+1)        + c_bar(irad))/TWO
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
                c_avg        = c_avg        + dr*(self%r(irad+1)*c_bar(irad+1)        + self%r(irad)*c_bar(irad))/TWO
                area = area + dr*(self%r(irad+1)+self%r(irad))/TWO
            end do

        end if

        density_avg  = density_avg  / area
        vel1_avg     = vel1_avg     / area
        vel2_avg     = vel2_avg     / area
        vel3_avg     = vel3_avg     / area
        pressure_avg = pressure_avg / area
        c_avg        = c_avg        / area

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
    subroutine analyze_bc_geometry(self,mesh,group_ID,bc_comm,side)
        class(nonlocal_nrbc_lindblad_base_t), intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh
        integer(ik),            intent(in)      :: group_ID
        type(mpi_comm),         intent(in)      :: bc_comm
        character(1),           intent(in)      :: side

        integer(ik) :: idomain_l, ielement_l, patch_ID, face_ID, iface, ierr, inode, local_nfaces_a, local_nfaces_b
        real(rk)    :: face_rmin, face_rmax, local_rmin, local_rmax, global_rmin, global_rmax, face_z,  &
                       face_thetamin, face_thetamax, local_thetamin, local_thetamax, global_thetamin, global_thetamax, &
                       local_zref_a, local_zref_b
        real(rk)    :: ref_nodes(8,3), physical_nodes(8,3)

        ! Zero face counter
        if (side=='A') local_nfaces_a = 0
        if (side=='B') local_nfaces_b = 0

        ! Search for min/max radius on local processor
        local_rmin     =  HUGE(1._rk) ! any radius will be smaller than this, so it is guarunteed to be reset.
        local_rmax     = -HUGE(1._rk) ! any radius will be larger than this, so it is guarunteed to be reset.
        local_thetamin =  HUGE(1._rk) ! any theta will be smaller than this, so it is guarunteed to be reset.
        local_thetamax = -HUGE(1._rk) ! any theta will be larger than this, so it is guarunteed to be reset.
        local_zref_a   =  HUGE(1._rk) ! default to huge negative number. We assume no actual coordinate will ever be this.
        local_zref_b   = -HUGE(1._rk) ! default to huge negative number. We assume no actual coordinate will ever be this.
        do patch_ID = 1,mesh%bc_patch_group(group_ID)%npatches()
            do face_ID = 1,mesh%bc_patch_group(group_ID)%patch(patch_ID)%nfaces()

                idomain_l  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
                ielement_l = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                iface      = mesh%bc_patch_group(group_ID)%patch(patch_ID)%iface(face_ID)


                ! Pick points on face to evaluate coordinate 
                ! extrema:
                !
                !   *---*---*
                !   |       |
                !   *       *
                !   |       |
                !   *---*---*
                !
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
                    call chidg_signal(FATAL,"outlet_giles_quasi3d_unsteady_HB: analyze_bc_geometry, invalid face index.")
                end if


                ! Handle geometry in two sets:
                !   'A' those faces with positive outward facing z-normal
                !   'B' those faces with negative outward facing z-normal
                if ( (side == 'A' .and. mesh%domain(idomain_l)%faces(ielement_l,iface)%unorm(1,3) > ZERO) .or. &
                     (side == 'B' .and. mesh%domain(idomain_l)%faces(ielement_l,iface)%unorm(1,3) < ZERO) )  then

                    ! Count number of faces in A/B sides
                    if (side=='A') local_nfaces_a = local_nfaces_a + 1
                    if (side=='B') local_nfaces_b = local_nfaces_b + 1

                    ! Evaluate physical coordinates on face edges
                    idomain_l  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%idomain_l()
                    ielement_l = mesh%bc_patch_group(group_ID)%patch(patch_ID)%ielement_l(face_ID)
                    do inode = 1,size(ref_nodes,1)
                        physical_nodes(inode,:) = mesh%domain(idomain_l)%elems(ielement_l)%physical_point(ref_nodes(inode,:),'Deformed')
                    end do !inode

                    ! Get face min/max radius
                    face_rmin     = minval(physical_nodes(:,1))
                    face_rmax     = maxval(physical_nodes(:,1))
                    face_thetamin = minval(physical_nodes(:,2))
                    face_thetamax = maxval(physical_nodes(:,2))

                    if (side=='A') face_z = minval(physical_nodes(:,3))
                    if (side=='B') face_z = maxval(physical_nodes(:,3))

                    ! Update processor-local value if new min/max values were found on the face
                    if (face_rmin < local_rmin)         local_rmin     = face_rmin
                    if (face_rmax > local_rmax)         local_rmax     = face_rmax
                    if (face_thetamin < local_thetamin) local_thetamin = face_thetamin
                    if (face_thetamax > local_thetamax) local_thetamax = face_thetamax

                    if (side=='A' .and. face_z < local_zref_a) local_zref_a = face_z
                    if (side=='B' .and. face_z > local_zref_b) local_zref_b = face_z

                end if

            end do !face_ID
        end do !patch_ID

        ! Reduce face counts
        call MPI_AllReduce(local_nfaces_a,self%nfaces_a,1,MPI_INTEGER4,MPI_SUM,bc_comm,ierr)
        call MPI_AllReduce(local_nfaces_b,self%nfaces_b,1,MPI_INTEGER4,MPI_SUM,bc_comm,ierr)

        ! Reduce processor local values to determine boundary-global min/max values
        call MPI_AllReduce(local_rmin,    global_rmin,    1,MPI_REAL8,MPI_MIN,bc_comm,ierr)
        call MPI_AllReduce(local_rmax,    global_rmax,    1,MPI_REAL8,MPI_MAX,bc_comm,ierr)
        call MPI_AllReduce(local_thetamin,global_thetamin,1,MPI_REAL8,MPI_MIN,bc_comm,ierr)
        call MPI_AllReduce(local_thetamax,global_thetamax,1,MPI_REAL8,MPI_MAX,bc_comm,ierr)
        
        ! Create radial stations
        if ( (side=='A' .and. self%nfaces_a > 0) .or. (side=='B' .and. self%nfaces_b > 0) ) then
            self%r = linspace(global_rmin + 0.0000001_rk,global_rmax-0.0000001_rk,self%nr)
        end if

        ! Compute z_ref
        if (side == 'A')then
            self%theta_ref_a = (global_thetamin + global_thetamax)/TWO
            if (self%nfaces_a /= 0) self%theta_ref = self%theta_ref_a
            call MPI_AllReduce(local_zref_a,self%z_ref_a,1,MPI_REAL8,MPI_MIN,bc_comm,ierr)
        else if (side == 'B') then
            self%theta_ref_b = (global_thetamin + global_thetamax)/TWO
            if (self%nfaces_b /= 0) self%theta_ref = self%theta_ref_b
            call MPI_AllReduce(local_zref_b,self%z_ref_b,1,MPI_REAL8,MPI_MAX,bc_comm,ierr)
        else
            call chidg_signal(FATAL,"bc_giles: analyze_bc_geometry invalid input for argument 'side'. 'A' or 'B'.")
        end if

        self%theta_ref = ZERO

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
    !!  NOTE: Assumes entire boundary is CONSTANT-Z
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
    subroutine initialize_fourier_discretization(self,mesh,group_ID,bc_comm,side)
        class(nonlocal_nrbc_lindblad_base_t), intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh
        integer(ik),            intent(in)      :: group_ID
        type(mpi_comm),         intent(in)      :: bc_comm
        character(1),           intent(in)      :: side

        integer(ik)                 :: nmodes, ncoeff, nradius, ntheta, idomain_l, ielement_l, iface, &
                                       iradius, itheta, ierr, noverset, pelem_ID, patch_ID, face_ID
        real(rk)                    :: dtheta, dtheta_n, midpoint(3), try_offset(3), node(3), z, donor_avg_z
        real(rk),       allocatable :: pitch(:), donor_z(:)
        character(:),   allocatable :: user_msg
        logical                     :: donor_found, found_face_on_side

        type(element_info_t),   allocatable :: donors(:,:)
        real(rk),               allocatable :: thetas(:,:)
        real(rk),               allocatable :: donor_nodes(:,:,:)
        real(rk)                            :: theta_ref

        type(rule_prefers_A)    :: donor_rule_A
        type(rule_prefers_B)    :: donor_rule_B
        class(multi_donor_rule_t), allocatable  :: multi_donor_rule

        ! If no faces on A/B, then the auxiliary grid for that side should not
        ! be constructed.
        if (side=='A' .and. self%nfaces_a==0) return
        if (side=='B' .and. self%nfaces_b==0) return


        ! Define Fourier discretization
        nmodes  = self%nfourier_space
        ncoeff  = 1 + (nmodes-1)*2
        nradius = size(self%r)
        ntheta  = ncoeff

        
        ! Initialize theta discretization parameters
        if (side == 'A') then
            pitch  = self%bcproperties%compute('Pitch A',time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])
            theta_ref = self%theta_ref_a
            z = self%z_ref_a
            allocate(multi_donor_rule,source=donor_rule_A,stat=ierr)
            if (ierr /= 0) call AllocationError
        else if (side == 'B') then
            pitch  = self%bcproperties%compute('Pitch B',time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])
            theta_ref = self%theta_ref_b
            z = self%z_ref_b
            allocate(multi_donor_rule,source=donor_rule_B,stat=ierr)
            if (ierr /= 0) call AllocationError
        end if
        dtheta = pitch(1)
        dtheta_n = dtheta/ntheta


        ! Construct theta discretization at each radius
        allocate(thetas(size(self%r),ntheta), donors(nradius,ntheta), donor_nodes(nradius,ntheta,3),  stat=ierr)
        if (ierr /= 0) call AllocationError
        do itheta = 1,ntheta
            thetas(:,itheta) = theta_ref + (itheta-1)*dtheta_n
        end do


        ! Donor search offset, if needed
        try_offset = [ZERO, -pitch(1), ZERO]

        ! For each radial station, initialized donor for each node in theta grid
        do iradius = 1,size(self%r)
            noverset = 0
            do itheta = 1,ntheta

                node = [self%r(iradius), thetas(iradius,itheta), z]

                ! Try processor-LOCAL elements
                call find_gq_donor(mesh,                                &
                                   node,                                &
                                   [ZERO,ZERO,ZERO],                    &
                                   face_info_constructor(0,0,0,0,0),    &   ! we don't really have a receiver face
                                   donors(iradius,itheta),              &
                                   donor_nodes(iradius,itheta,1:3),     &
                                   donor_found,                         &
                                   multi_donor_rule=multi_donor_rule)


                ! Reject if not correct face_normal for side: do this by checking the element center
                ! against the z-constant value.
                if (donor_found) then
                    donor_z = mesh%domain(donors(iradius,itheta)%idomain_l)%elems(donors(iradius,itheta)%ielement_l)%node_coords_def(:,3)
                    donor_avg_z = sum(donor_z)/size(donor_z)
                    if ( (side == 'A') .and. (donor_avg_z > z) ) donor_found = .false.
                    if ( (side == 'B') .and. (donor_avg_z < z) ) donor_found = .false.
                end if


                ! Try LOCAL elements with try_offset if still not found 
                if ( .not. donor_found ) then
                    call find_gq_donor(mesh,                                &
                                       node,                                &
                                       try_offset,                          &
                                       face_info_constructor(0,0,0,0,0),    &   ! we don't really have a receiver face
                                       donors(iradius,itheta),              &
                                       donor_nodes(iradius,itheta,1:3),     &
                                       donor_found,                         &
                                       multi_donor_rule=multi_donor_rule)


                    ! Reject if not correct face_normal for side: do this by checking the element center
                    ! against the z-constant value.
                    if (donor_found) then
                        donor_z = mesh%domain(donors(iradius,itheta)%idomain_l)%elems(donors(iradius,itheta)%ielement_l)%node_coords_def(:,3)
                        donor_avg_z = sum(donor_z)/size(donor_z)
                        if ( (side == 'A') .and. (donor_avg_z > z) ) donor_found = .false.
                        if ( (side == 'B') .and. (donor_avg_z < z) ) donor_found = .false.
                    end if

                    if (donor_found) then
                        noverset=noverset+1
                        thetas(iradius,itheta) = thetas(iradius,itheta) - pitch(1)
                    end if

                end if


                ! Try PARALLEL_ELEMENTS if donor not found amongst local elements
                if (.not. donor_found) then
                    call find_gq_donor_parallel(mesh,                               &
                                                node,                               &
                                                [ZERO,ZERO,ZERO],                   &
                                                face_info_constructor(0,0,0,0,0),   &   ! we don't really have a receiver face
                                                donors(iradius,itheta),             &
                                                donor_nodes(iradius,itheta,1:3),    &
                                                donor_found)

                    ! Reject if not correct face_normal for side: do this by checking the element center
                    ! against the z-constant value.
                    if (donor_found) then
                        pelem_ID = donors(iradius,itheta)%pelem_ID
                        donor_z = mesh%parallel_element(pelem_ID)%node_coords_def(:,3)
                        donor_avg_z = sum(donor_z)/size(donor_z)
                        if ( (side == 'A') .and. (donor_avg_z > z) ) donor_found = .false.
                        if ( (side == 'B') .and. (donor_avg_z < z) ) donor_found = .false.
                    end if
                end if

                
                ! Try PARALLEL_ELEMENTS with try_offset if still not found 
                if ( .not. donor_found ) then
                    call find_gq_donor_parallel(mesh,                               &
                                                node,                               &
                                                try_offset,                         &
                                                face_info_constructor(0,0,0,0,0),   &   ! we don't really have a receiver face
                                                donors(iradius,itheta),             &
                                                donor_nodes(iradius,itheta,1:3),    &
                                                donor_found)

                    ! Reject if not correct face_normal for side: do this by checking the element center
                    ! against the z-constant value.
                    if (donor_found) then
                        pelem_ID = donors(iradius,itheta)%pelem_ID
                        donor_z = mesh%parallel_element(pelem_ID)%node_coords_def(:,3)
                        donor_avg_z = sum(donor_z)/size(donor_z)
                        if ( (side == 'A') .and. (donor_avg_z > z) ) donor_found = .false.
                        if ( (side == 'B') .and. (donor_avg_z < z) ) donor_found = .false.
                    end if

                    if (donor_found) then
                        noverset=noverset+1
                        thetas(iradius,itheta) = thetas(iradius,itheta) - pitch(1)
                    end if
                end if 


                ! Abort if we didn't find a donor
                user_msg = "bc_nonlocal_nrbc_lindblad_base%initialize_fourier_discretization: &
                            no donor element found for Fourier discretization node."
                if (.not. donor_found) call chidg_signal(FATAL,user_msg)

            end do !itheta

            ! Shift arrays so that we start with the theta_min point
            donors(iradius,:)        = cshift(donors(iradius,:),        -noverset, dim=1)
            donor_nodes(iradius,:,:) = cshift(donor_nodes(iradius,:,:), -noverset, dim=1)
            thetas(iradius,:)        = cshift(thetas(iradius,:),        -noverset, dim=1)

        end do !iradius


        ! Store to correct side
        if (side == 'A') then
            self%donor_a      = donors
            self%donor_node_a = donor_nodes
            self%theta_a      = thetas
        else if (side == 'B') then
            self%donor_b      = donors
            self%donor_node_b = donor_nodes
            self%theta_b      = thetas
        end if


        ! If other side doesn't have any faces, copy current side discretiation
        if (side=='A' .and. self%nfaces_b==0) then
            self%theta_b = self%theta_a
        else if (side=='B' .and. self%nfaces_a==0) then
            self%theta_a = self%theta_b
        end if

    end subroutine initialize_fourier_discretization
    !**************************************************************************************





    !>  Return the correct temporal frequency for the time index.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/28/2018
    !!
    !--------------------------------------------------------------------------------------
    function get_omega(worker,itime) result(omega)
        type(chidg_worker_t),   intent(in)  :: worker
        integer(ik),            intent(in)  :: itime

        real(rk) :: omega
        integer(ik) :: ntime

        ntime = worker%time_manager%ntime

        if (itime == 1) then
            omega = ZERO
        else if (itime <= ((ntime-1)/2 + 1)) then
            omega = worker%time_manager%freqs(itime-1)
        else
            omega = -worker%time_manager%freqs(ntime-itime+1)
        end if



!        if (itime == 1) then
!            omega = ZERO
!        else if (itime == 2) then
!            omega = worker%time_manager%freqs(1)
!        else if (itime == 3) then
!            omega = -worker%time_manager%freqs(1)
!        end if


    end function get_omega
    !**************************************************************************************


    !>  Return the correct spatial frequency for the space index.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/28/2018
    !!
    !--------------------------------------------------------------------------------------
    function get_lm(itheta,ntheta) result(lm)
        integer(ik),    intent(in)  :: itheta
        integer(ik),    intent(in)  :: ntheta

        real(rk) :: lm

        if (itheta <= ((ntheta-1)/2 + 1)) then
            lm = real(itheta-1,rk)  ! positive frequencies
        else
            lm = -real(ntheta-itheta+1,rk) ! negative frequencies
        end if

    end function get_lm
    !**************************************************************************************
    

    !>  Return the side of the interface that the element currently pointed to by the
    !!  chidg_worker belongs to. 'A' or 'B'
    !!
    !!  Determined by the direction of the third component of the unit normal for 
    !!  the current face.
    !!
    !!  Returns: side = 'A' or 'B'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/13/2018
    !!
    !--------------------------------------------------------------------------------------
    function get_face_side(self,worker) result(side)
        class(nonlocal_nrbc_lindblad_base_t), intent(in)  :: self
        type(chidg_worker_t),   intent(in)  :: worker

        real(rk), allocatable, dimension(:) :: unorm3
        character(1)    :: side

        unorm3 = worker%unit_normal(3)

        if (unorm3(1) > ZERO) side = 'A'
        if (unorm3(1) < ZERO) side = 'B'

    end function get_face_side
    !**************************************************************************************

    !>  When multiple donors exist during construction of the Fourier interface 
    !!  interpolations, for side 'A', prefer donors with the smaller z-coordinate value
    !!
    !--------------------------------------------------------------------------------------
    function select_donor_A(mesh,donors,candidate_domains_g,candidate_domains_l,candidate_elements_g,candidate_elements_l) result(donor_index)
        type(mesh_t),       intent(in)  :: mesh
        type(ivector_t),    intent(in)  :: donors
        type(ivector_t),    intent(in)  :: candidate_domains_g
        type(ivector_t),    intent(in)  :: candidate_domains_l
        type(ivector_t),    intent(in)  :: candidate_elements_g
        type(ivector_t),    intent(in)  :: candidate_elements_l
        integer(ik) :: donor_index

        real(rk), allocatable   :: donor_z(:)
        integer(ik) :: idonor

        ! Get index of domain with minimum volume
        allocate(donor_z(donors%size()))
        do idonor = 1,donors%size()
            donor_z(idonor) = minval(mesh%domain(candidate_domains_l%at(donors%at(idonor)))%elems(candidate_elements_l%at(donors%at(idonor)))%node_coords_def(:,3))
        end do 
        donor_index = minloc(donor_z,1)

    end function select_donor_A

    !>  When multiple donors exist during construction of the Fourier interface 
    !!  interpolations, for side 'B', prefer donors with the larger z-coordinate value
    !!
    !--------------------------------------------------------------------------------------
    function select_donor_B(mesh,donors,candidate_domains_g,candidate_domains_l,candidate_elements_g,candidate_elements_l) result(donor_index)
        type(mesh_t),       intent(in)  :: mesh
        type(ivector_t),    intent(in)  :: donors
        type(ivector_t),    intent(in)  :: candidate_domains_g
        type(ivector_t),    intent(in)  :: candidate_domains_l
        type(ivector_t),    intent(in)  :: candidate_elements_g
        type(ivector_t),    intent(in)  :: candidate_elements_l
        integer(ik) :: donor_index

        real(rk), allocatable   :: donor_z(:)
        integer(ik) :: idonor

        ! Get index of domain with minimum volume
        allocate(donor_z(donors%size()))
        do idonor = 1,donors%size()
            donor_z(idonor) = maxval(mesh%domain(candidate_domains_l%at(donors%at(idonor)))%elems(candidate_elements_l%at(donors%at(idonor)))%node_coords_def(:,3))
        end do 
        donor_index = maxloc(donor_z,1)

    end function select_donor_B




end module bc_nonlocal_nrbc_lindblad_base
