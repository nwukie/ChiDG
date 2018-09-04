module bc_state_outlet_characteristic_quasi3d_unsteady_HB
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF, ME, CYLINDRICAL,    &
                                      XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, PI, NO_PROC
    use mod_fluid,              only: gam
    use mod_interpolation,      only: interpolate_linear, interpolate_linear_ad
    use mod_gridspace,          only: linspace
    use mod_dft,                only: dft, idft_eval
    use mod_chimera,            only: find_gq_donor, find_gq_donor_parallel

    use type_point,             only: point_t
    use type_mesh,              only: mesh_t
    use type_bc_state,          only: bc_state_t
    use bc_state_fluid_averaging,   only: bc_fluid_averaging_t
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
    type, public, extends(bc_fluid_averaging_t) :: outlet_characteristic_quasi3d_unsteady_HB_t

        integer(ik) :: nr = 10
        integer(ik) :: nfourier_space = 14

        real(rk),   allocatable :: r(:)
        real(rk),   allocatable :: theta(:,:)   ! (nr,ntheta)
        real(rk)                :: theta_ref    

        type(element_info_t),   allocatable :: donor(:,:)
        real(rk),               allocatable :: donor_node(:,:,:)

    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: init_bc_postcomm     ! Implement specialized initialization
        procedure   :: compute_bc_state     ! boundary condition function implementation

        procedure   :: compute_temporal_dft
        procedure   :: compute_spatial_dft
        procedure   :: apply_nonreflecting_condition
        procedure   :: compute_steady_nrbc
        procedure   :: compute_steady_decomposition
        procedure   :: compute_boundary_average
        procedure   :: analyze_bc_geometry
        procedure   :: initialize_fourier_discretization

    end type outlet_characteristic_quasi3d_unsteady_HB_t
    !*********************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. average_pressure 
    !!  @date   2/8/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(outlet_characteristic_quasi3d_unsteady_HB_t),   intent(inout) :: self

        ! Set name, family
        call self%set_name('Outlet - Characteristic Quasi3D Unsteady HB')
        call self%set_family('Outlet')

        ! Add functions
        call self%bcproperties%add('Average Pressure',    'Required')
        call self%bcproperties%add('Pitch',               'Required')
        call self%bcproperties%add('Spatial Periodicity', 'Required')

    end subroutine init
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
        class(outlet_characteristic_quasi3d_unsteady_HB_t),   intent(inout)   :: self
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
    !!  @date   2/8/2018
    !!
    !!  @param[in]      worker  Interface for geometry, cache, integration, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_comm)
        class(outlet_characteristic_quasi3d_unsteady_HB_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop
        type(mpi_comm),                             intent(in)      :: bc_comm


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                                                      &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc, pressure_bc, vel1_bc, vel2_bc, vel3_bc,   &
            density_bc_tmp, vel1_bc_tmp, vel2_bc_tmp, vel3_bc_tmp, pressure_bc_tmp,                     &
            grad1_density_m, grad1_mom1_m, grad1_mom2_m, grad1_mom3_m, grad1_energy_m,                  &
            grad2_density_m, grad2_mom1_m, grad2_mom2_m, grad2_mom3_m, grad2_energy_m,                  &
            grad3_density_m, grad3_mom1_m, grad3_mom2_m, grad3_mom3_m, grad3_energy_m, expect_zero
!            c1,    c2,    c3,    c4,    c5,                                                             &
!            c1_3d, c2_3d, c3_3d, c4_3d, c5_3d,                                                          &
!            c1_1d, c2_1d, c3_1d, c4_1d, c5_1d,                                                          &
!            density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, c_bar,                             &
!            ddensity, dvel1, dvel2, dvel3, dpressure, expect_zero

        type(AD_D), allocatable, dimension(:,:) ::                                  &
            density_t_real, vel1_t_real, vel2_t_real, vel3_t_real, pressure_t_real, &
            density_t_imag, vel1_t_imag, vel2_t_imag, vel3_t_imag, pressure_t_imag


        type(AD_D), allocatable, dimension(:,:,:) ::                                            &
            density_Ft_real, vel1_Ft_real, vel2_Ft_real, vel3_Ft_real, pressure_Ft_real,        &
            density_Ft_imag, vel1_Ft_imag, vel2_Ft_imag, vel3_Ft_imag, pressure_Ft_imag,        &
            density_Fts_real, vel1_Fts_real, vel2_Fts_real, vel3_Fts_real, pressure_Fts_real,   &
            density_Fts_imag, vel1_Fts_imag, vel2_Fts_imag, vel3_Fts_imag, pressure_Fts_imag,   &
            density_Fts_real_gq, vel1_Fts_real_gq, vel2_Fts_real_gq, vel3_Fts_real_gq, pressure_Fts_real_gq,   &
            density_Fts_imag_gq, vel1_Fts_imag_gq, vel2_Fts_imag_gq, vel3_Fts_imag_gq, pressure_Fts_imag_gq


        real(rk),       allocatable, dimension(:)   :: p_user, r, pitch
        real(rk)                                    :: theta_offset
        type(point_t),  allocatable                 :: coords(:)
        integer                                     :: i, ngq, ivec, imode, itheta, itime, iradius, nmodes, ierr, igq


        ! Get back pressure from function.
        p_user = self%bcproperties%compute('Average Pressure',worker%time(),worker%coords())
        pitch  = self%bcproperties%compute('Pitch',           worker%time(),worker%coords())

!        ! Interpolate interior solution to face quadrature nodes
!        density_m = worker%get_field('Density'   , 'value', 'face interior')
!        mom1_m    = worker%get_field('Momentum-1', 'value', 'face interior')
!        mom2_m    = worker%get_field('Momentum-2', 'value', 'face interior')
!        mom3_m    = worker%get_field('Momentum-3', 'value', 'face interior')
!        energy_m  = worker%get_field('Energy'    , 'value', 'face interior')

        grad1_density_m = worker%get_field('Density'   , 'grad1', 'face interior')
        grad2_density_m = worker%get_field('Density'   , 'grad2', 'face interior')
        grad3_density_m = worker%get_field('Density'   , 'grad3', 'face interior')

        grad1_mom1_m    = worker%get_field('Momentum-1', 'grad1', 'face interior')
        grad2_mom1_m    = worker%get_field('Momentum-1', 'grad2', 'face interior')
        grad3_mom1_m    = worker%get_field('Momentum-1', 'grad3', 'face interior')

        grad1_mom2_m    = worker%get_field('Momentum-2', 'grad1', 'face interior')
        grad2_mom2_m    = worker%get_field('Momentum-2', 'grad2', 'face interior')
        grad3_mom2_m    = worker%get_field('Momentum-2', 'grad3', 'face interior')

        grad1_mom3_m    = worker%get_field('Momentum-3', 'grad1', 'face interior')
        grad2_mom3_m    = worker%get_field('Momentum-3', 'grad2', 'face interior')
        grad3_mom3_m    = worker%get_field('Momentum-3', 'grad3', 'face interior')
        
        grad1_energy_m  = worker%get_field('Energy'    , 'grad1', 'face interior')
        grad2_energy_m  = worker%get_field('Energy'    , 'grad2', 'face interior')
        grad3_energy_m  = worker%get_field('Energy'    , 'grad3', 'face interior')


        ! Store boundary gradient state. Grad(Q_bc). Do this here, before we
        ! compute any transformations for cylindrical.
        call worker%store_bc_state('Density'   , grad1_density_m, 'grad1')
        call worker%store_bc_state('Density'   , grad2_density_m, 'grad2')
        call worker%store_bc_state('Density'   , grad3_density_m, 'grad3')
                                                
        call worker%store_bc_state('Momentum-1', grad1_mom1_m,    'grad1')
        call worker%store_bc_state('Momentum-1', grad2_mom1_m,    'grad2')
        call worker%store_bc_state('Momentum-1', grad3_mom1_m,    'grad3')
                                                
        call worker%store_bc_state('Momentum-2', grad1_mom2_m,    'grad1')
        call worker%store_bc_state('Momentum-2', grad2_mom2_m,    'grad2')
        call worker%store_bc_state('Momentum-2', grad3_mom2_m,    'grad3')
                                                
        call worker%store_bc_state('Momentum-3', grad1_mom3_m,    'grad1')
        call worker%store_bc_state('Momentum-3', grad2_mom3_m,    'grad2')
        call worker%store_bc_state('Momentum-3', grad3_mom3_m,    'grad3')
                                                
        call worker%store_bc_state('Energy'    , grad1_energy_m,  'grad1')
        call worker%store_bc_state('Energy'    , grad2_energy_m,  'grad2')
        call worker%store_bc_state('Energy'    , grad3_energy_m,  'grad3')


!        ! Account for cylindrical. Get tangential momentum from angular momentum.
!        r = worker%coordinate('1','boundary')
!        if (worker%coordinate_system() == 'Cylindrical') then
!            mom2_m = mom2_m / r
!            grad1_mom2_m = (grad1_mom2_m/r) - mom2_m/r
!            grad2_mom2_m = (grad2_mom2_m/r)
!            grad3_mom2_m = (grad3_mom2_m/r)
!        end if
!
!        ! Compute velocity and pressure
!        vel1_m = mom1_m/density_m
!        vel2_m = mom2_m/density_m
!        vel3_m = mom3_m/density_m
!        pressure_m = worker%get_field('Pressure', 'value', 'face interior')


        ! Compute Fourier decomposition of temporal data at points
        ! on the spatial transform grid.
        !   : U_Ft(nradius,ntheta,ntime)
        call self%compute_temporal_dft(worker,bc_comm,                        &
                                       density_Ft_real,  density_Ft_imag,     &
                                       vel1_Ft_real,     vel1_Ft_imag,        &
                                       vel2_Ft_real,     vel2_Ft_imag,        &
                                       vel3_Ft_real,     vel3_Ft_imag,        &
                                       pressure_Ft_real, pressure_Ft_imag)


        ! Compute Fourier decomposition in theta at set of radial 
        ! stations for each temporal mode:
        !   : U_Fts(nradius,ntheta,ntime)
        call self%compute_spatial_dft(worker,bc_comm,                         &
                                      density_Ft_real,   density_Ft_imag,     &
                                      vel1_Ft_real,      vel1_Ft_imag,        &
                                      vel2_Ft_real,      vel2_Ft_imag,        &
                                      vel3_Ft_real,      vel3_Ft_imag,        &
                                      pressure_Ft_real,  pressure_Ft_imag,    &
                                      density_Fts_real,  density_Fts_imag,    &
                                      vel1_Fts_real,     vel1_Fts_imag,       &
                                      vel2_Fts_real,     vel2_Fts_imag,       &
                                      vel3_Fts_real,     vel3_Fts_imag,       &
                                      pressure_Fts_real, pressure_Fts_imag)


        !! Get spatio-temporal average at each radial station
        !density_bar  = density_Fts_real(:,1,1)
        !vel1_bar     = vel1_Fts_real(:,1,1)
        !vel2_bar     = vel2_Fts_real(:,1,1)
        !vel3_bar     = vel3_Fts_real(:,1,1)
        !pressure_bar = pressure_Fts_real(:,1,1)

        call compute_steady_nrbc(self,worker,bc_comm,       &
                                 density_Fts_real(:,:,1),   &
                                 density_Fts_imag(:,:,1),   &
                                 vel1_Fts_real(:,:,1),      &
                                 vel1_Fts_imag(:,:,1),      &
                                 vel2_Fts_real(:,:,1),      &
                                 vel2_Fts_imag(:,:,1),      &
                                 vel3_Fts_real(:,:,1),      &
                                 vel3_Fts_imag(:,:,1),      &
                                 pressure_Fts_real(:,:,1),  &
                                 pressure_Fts_imag(:,:,1))


!        !
!        ! Apply nonreflecting analysis to the spatio-temporal Fourier decomposition
!        ! at each radial station.
!        !
!        call self%apply_nonreflecting_condition(worker,bc_comm,                         &
!                                                density_Fts_real,  density_Fts_imag,    &
!                                                vel1_Fts_real,     vel1_Fts_imag,       &
!                                                vel2_Fts_real,     vel2_Fts_imag,       &
!                                                vel3_Fts_real,     vel3_Fts_imag,       &
!                                                pressure_Fts_real, pressure_Fts_imag)



        ! Interpolate spatio-temporal Fourier coefficients to quadrature nodes
        ! linear interpolation between radial coordinates.
        !allocate(density_Fts_real(nradius,ntheta,ntime), density_Fts_imag(nradius,ntheta,ntime),   &
        coords = worker%coords()
        allocate(density_Fts_real_gq( size(coords),size(density_Fts_real,2),size(density_Fts_real,3)), &
                 density_Fts_imag_gq( size(coords),size(density_Fts_real,2),size(density_Fts_real,3)), &
                 vel1_Fts_real_gq(    size(coords),size(density_Fts_real,2),size(density_Fts_real,3)), &
                 vel1_Fts_imag_gq(    size(coords),size(density_Fts_real,2),size(density_Fts_real,3)), &
                 vel2_Fts_real_gq(    size(coords),size(density_Fts_real,2),size(density_Fts_real,3)), &
                 vel2_Fts_imag_gq(    size(coords),size(density_Fts_real,2),size(density_Fts_real,3)), &
                 vel3_Fts_real_gq(    size(coords),size(density_Fts_real,2),size(density_Fts_real,3)), &
                 vel3_Fts_imag_gq(    size(coords),size(density_Fts_real,2),size(density_Fts_real,3)), &
                 pressure_Fts_real_gq(size(coords),size(density_Fts_real,2),size(density_Fts_real,3)), &
                 pressure_Fts_imag_gq(size(coords),size(density_Fts_real,2),size(density_Fts_real,3)), &
                 stat=ierr)
        if (ierr /= 0) call AllocationError
        density_Fts_real_gq  = ZERO*density_Fts_real(1,1,1)
        density_Fts_imag_gq  = ZERO*density_Fts_real(1,1,1)
        vel1_Fts_real_gq     = ZERO*density_Fts_real(1,1,1)
        vel1_Fts_imag_gq     = ZERO*density_Fts_real(1,1,1)
        vel2_Fts_real_gq     = ZERO*density_Fts_real(1,1,1)
        vel2_Fts_imag_gq     = ZERO*density_Fts_real(1,1,1)
        vel3_Fts_real_gq     = ZERO*density_Fts_real(1,1,1)
        vel3_Fts_imag_gq     = ZERO*density_Fts_real(1,1,1)
        pressure_Fts_real_gq = ZERO*density_Fts_real(1,1,1)
        pressure_Fts_imag_gq = ZERO*density_Fts_real(1,1,1)
        do igq = 1,size(coords)
            do itheta = 1,size(density_Fts_real,2)
                do itime = 1,size(density_Fts_real,3)
                    density_Fts_real_gq( igq,itheta,itime) = interpolate_linear_ad(self%r,density_Fts_real(:,itheta,itime), coords(igq)%c1_)
                    density_Fts_imag_gq( igq,itheta,itime) = interpolate_linear_ad(self%r,density_Fts_imag(:,itheta,itime), coords(igq)%c1_)
                    vel1_Fts_real_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel1_Fts_real(:,itheta,itime),    coords(igq)%c1_)
                    vel1_Fts_imag_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel1_Fts_imag(:,itheta,itime),    coords(igq)%c1_)
                    vel2_Fts_real_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel2_Fts_real(:,itheta,itime),    coords(igq)%c1_)
                    vel2_Fts_imag_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel2_Fts_imag(:,itheta,itime),    coords(igq)%c1_)
                    vel3_Fts_real_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel3_Fts_real(:,itheta,itime),    coords(igq)%c1_)
                    vel3_Fts_imag_gq(    igq,itheta,itime) = interpolate_linear_ad(self%r,vel3_Fts_imag(:,itheta,itime),    coords(igq)%c1_)
                    pressure_Fts_real_gq(igq,itheta,itime) = interpolate_linear_ad(self%r,pressure_Fts_real(:,itheta,itime),coords(igq)%c1_)
                    pressure_Fts_imag_gq(igq,itheta,itime) = interpolate_linear_ad(self%r,pressure_Fts_imag(:,itheta,itime),coords(igq)%c1_)
                end do !itime
            end do !itheta
        end do !igq

!        ! Zero out Fourier modes for Zero-th temporal mode. We will handle the time-average
!        ! with the standard Steady Giles approach
!        density_Fts_real_gq(:,:,1)  = ZERO
!        density_Fts_imag_gq(:,:,1)  = ZERO
!        vel1_Fts_real_gq(:,:,1)     = ZERO
!        vel1_Fts_imag_gq(:,:,1)     = ZERO
!        vel2_Fts_real_gq(:,:,1)     = ZERO
!        vel2_Fts_imag_gq(:,:,1)     = ZERO
!        vel3_Fts_real_gq(:,:,1)     = ZERO
!        vel3_Fts_imag_gq(:,:,1)     = ZERO
!        pressure_Fts_real_gq(:,:,1) = ZERO
!        pressure_Fts_imag_gq(:,:,1) = ZERO
!
!        ! Impose constant pressure
!        pressure_Fts_real_gq(:,:,:) = ZERO
!        pressure_Fts_imag_gq(:,:,:) = ZERO
!        pressure_Fts_real_gq(:,1,1) = 100000._rk

        ! Reconstruct primitive variables at quadrature nodes from absorbing Fourier modes
        ! via inverse transform.

        ! Inverse DFT of spatio-temporal Fourier modes to give temporal Fourier
        ! modes at quadrature nodes.
        expect_zero = [AD_D(1)]
        density_t_real  = ZERO*density_Fts_real_gq(:,1,:)
        vel1_t_real     = ZERO*density_Fts_real_gq(:,1,:)
        vel2_t_real     = ZERO*density_Fts_real_gq(:,1,:)
        vel3_t_real     = ZERO*density_Fts_real_gq(:,1,:)
        pressure_t_real = ZERO*density_Fts_real_gq(:,1,:)
        density_t_imag  = ZERO*density_Fts_real_gq(:,1,:)
        vel1_t_imag     = ZERO*density_Fts_real_gq(:,1,:)
        vel2_t_imag     = ZERO*density_Fts_real_gq(:,1,:)
        vel3_t_imag     = ZERO*density_Fts_real_gq(:,1,:)
        pressure_t_imag = ZERO*density_Fts_real_gq(:,1,:)
        do igq = 1,size(coords)
            do itime = 1,size(density_Fts_real_gq,3)
                theta_offset = coords(igq)%c2_ - self%theta_ref
                ! **** WARNING: probably want ipdft_eval here ****
                call idft_eval(density_Fts_real_gq(igq,:,itime),    &
                               density_Fts_imag_gq(igq,:,itime),    &
                               [theta_offset]/pitch(1),             &
                               density_t_real(igq:igq,itime),       &
                               density_t_imag(igq:igq,itime))

                call idft_eval(vel1_Fts_real_gq(igq,:,itime),   &
                               vel1_Fts_imag_gq(igq,:,itime),   &
                               [theta_offset]/pitch(1),         &
                               vel1_t_real(igq:igq,itime),      &
                               vel1_t_imag(igq:igq,itime))

                call idft_eval(vel2_Fts_real_gq(igq,:,itime),   &
                               vel2_Fts_imag_gq(igq,:,itime),   &
                               [theta_offset]/pitch(1),         &
                               vel2_t_real(igq:igq,itime),      &
                               vel2_t_imag(igq:igq,itime))

                call idft_eval(vel3_Fts_real_gq(igq,:,itime),   &
                               vel3_Fts_imag_gq(igq,:,itime),   &
                               [theta_offset]/pitch(1),         &
                               vel3_t_real(igq:igq,itime),      &
                               vel3_t_imag(igq:igq,itime))

                call idft_eval(pressure_Fts_real_gq(igq,:,itime),   &
                               pressure_Fts_imag_gq(igq,:,itime),   &
                               [theta_offset]/pitch(1),             &
                               pressure_t_real(igq:igq,itime),      &
                               pressure_t_imag(igq:igq,itime))
            end do !itime
        end do !igq

        density_bc  = ZERO*density_Fts_real_gq(:,1,1)
        vel1_bc     = ZERO*density_Fts_real_gq(:,1,1)
        vel2_bc     = ZERO*density_Fts_real_gq(:,1,1)
        vel3_bc     = ZERO*density_Fts_real_gq(:,1,1)
        pressure_bc = ZERO*density_Fts_real_gq(:,1,1)

        ! Inverse DFT of temporal Fourier modes to give primitive variables
        ! at quarature nodes for the current time instance.
        density_bc_tmp  = [ZERO*density_Fts_real_gq(1,1,1)]
        vel1_bc_tmp     = [ZERO*density_Fts_real_gq(1,1,1)]
        vel2_bc_tmp     = [ZERO*density_Fts_real_gq(1,1,1)]
        vel3_bc_tmp     = [ZERO*density_Fts_real_gq(1,1,1)]
        pressure_bc_tmp = [ZERO*density_Fts_real_gq(1,1,1)]
        do igq = 1,size(coords)
            ! **** WARNING: probably want ipdft_eval here ****
            call idft_eval(density_t_real(igq,:),   &
                           density_t_imag(igq,:),   &
                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
                           density_bc_tmp,          &
                           expect_zero)

            call idft_eval(vel1_t_real(igq,:),      &
                           vel1_t_imag(igq,:),      &
                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
                           vel1_bc_tmp,        &
                           expect_zero)

            call idft_eval(vel2_t_real(igq,:),      &
                           vel2_t_imag(igq,:),      &
                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
                           vel2_bc_tmp,        &
                           expect_zero)

            call idft_eval(vel3_t_real(igq,:),      &
                           vel3_t_imag(igq,:),      &
                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
                           vel3_bc_tmp,        &
                           expect_zero)

            call idft_eval(pressure_t_real(igq,:),  &
                           pressure_t_imag(igq,:),  &
                           [real(worker%itime-1,rk)/real(worker%time_manager%ntime,rk)],    &
                           pressure_bc_tmp,    &
                           expect_zero)

            ! Accumulate contribution from unsteady modes
            density_bc(igq)  = density_bc(igq)  + density_bc_tmp(1)
            vel1_bc(igq)     = vel1_bc(igq)     + vel1_bc_tmp(1)
            vel2_bc(igq)     = vel2_bc(igq)     + vel2_bc_tmp(1)
            vel3_bc(igq)     = vel3_bc(igq)     + vel3_bc_tmp(1)
            pressure_bc(igq) = pressure_bc(igq) + pressure_bc_tmp(1)
        end do



        !
        ! Form conserved variables
        !
        density_bc = density_bc
        mom1_bc    = density_bc*vel1_bc
        mom2_bc    = density_bc*vel2_bc
        mom3_bc    = density_bc*vel3_bc
        energy_bc  = pressure_bc/(gam - ONE)  + (density_bc*HALF)*(vel1_bc*vel1_bc + vel2_bc*vel2_bc + vel3_bc*vel3_bc)


        !
        ! Account for cylindrical. Convert tangential momentum back to angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_bc = mom2_bc * r
        end if


        !
        ! Store boundary condition state. Q_bc
        !
        call worker%store_bc_state('Density'   , density_bc, 'value')
        call worker%store_bc_state('Momentum-1', mom1_bc,    'value')
        call worker%store_bc_state('Momentum-2', mom2_bc,    'value')
        call worker%store_bc_state('Momentum-3', mom3_bc,    'value')
        call worker%store_bc_state('Energy'    , energy_bc,  'value')



    end subroutine compute_bc_state
    !*********************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/12/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_temporal_dft(self,worker,bc_comm,                    &
                                    density_Ft_real,  density_Ft_imag,      &
                                    vel1_Ft_real,     vel1_Ft_imag,         &
                                    vel2_Ft_real,     vel2_Ft_imag,         &
                                    vel3_Ft_real,     vel3_Ft_imag,         &
                                    pressure_Ft_real, pressure_Ft_imag)
        class(outlet_characteristic_quasi3d_unsteady_HB_t),  intent(inout)   :: self
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


        type(AD_D), allocatable,    dimension(:)    ::                                          &
            density_real_tmp, vel1_real_tmp, vel2_real_tmp, vel3_real_tmp, pressure_real_tmp,   &
            density_imag_tmp, vel1_imag_tmp, vel2_imag_tmp, vel3_imag_tmp, pressure_imag_tmp


        type(AD_D), allocatable, dimension(:,:)     ::  &
            density, mom1, mom2, mom3, energy, vel1, vel2, vel3, pressure

        type(AD_D)  :: density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, c_bar

        integer(ik) :: nradius, ntheta, iradius, itheta, imode, itime, ntime, ierr



        ! Define Fourier space discretization to determine
        ! number of theta-samples are being taken
        ntheta  = size(self%theta,2)
        nradius = size(self%r)
        ntime   = worker%mesh%domain(worker%element_info%idomain_l)%elems(worker%element_info%ielement_l)%ntime

        ! Allocate storage for discrete time instances
        allocate(density( ntheta,ntime),  &
                 mom1(    ntheta,ntime),  &
                 mom2(    ntheta,ntime),  &
                 mom3(    ntheta,ntime),  &
                 vel1(    ntheta,ntime),  &
                 vel2(    ntheta,ntime),  &
                 vel3(    ntheta,ntime),  &
                 energy(  ntheta,ntime),  &
                 pressure(ntheta,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError
        
        ! Allocate storage for temporal dft
        allocate(density_Ft_real( nradius,ntheta,ntime), density_Ft_imag( nradius,ntheta,ntime),  &
                 vel1_Ft_real(    nradius,ntheta,ntime), vel1_Ft_imag(    nradius,ntheta,ntime),  &
                 vel2_Ft_real(    nradius,ntheta,ntime), vel2_Ft_imag(    nradius,ntheta,ntime),  &
                 vel3_Ft_real(    nradius,ntheta,ntime), vel3_Ft_imag(    nradius,ntheta,ntime),  &
                 pressure_Ft_real(nradius,ntheta,ntime), pressure_Ft_imag(nradius,ntheta,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Perform Fourier decomposition at each radial station.
        do iradius = 1,nradius

            ! Construct theta discretization
            do itime = 1,ntime
                ! Interpolate solution to physical_nodes at current radial station: [ntheta]
                density(:,itime) = worker%interpolate_field_general('Density',    donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)
                mom1(:,itime)    = worker%interpolate_field_general('Momentum-1', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)
                mom2(:,itime)    = worker%interpolate_field_general('Momentum-2', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)
                mom3(:,itime)    = worker%interpolate_field_general('Momentum-3', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)
                energy(:,itime)  = worker%interpolate_field_general('Energy',     donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:), itime=itime)

                if (worker%coordinate_system() == 'Cylindrical') then
                    mom2(:,itime) = mom2(:,itime)/self%r(iradius)  ! convert to tangential momentum
                end if
            end do

            ! Compute velocities and pressure at each time
            vel1 = mom1/density
            vel2 = mom2/density
            vel3 = mom3/density
            pressure = (gam-ONE)*(energy - HALF*( mom1*mom1 + mom2*mom2 + mom3*mom3 )/density )


            ! Temporal dft at each theta location
            do itheta = 1,ntheta
                call dft(density(itheta,:),  ZERO*density(itheta,:),  density_real_tmp,  density_imag_tmp )
                call dft(vel1(itheta,:),     ZERO*vel1(itheta,:),     vel1_real_tmp,     vel1_imag_tmp    )
                call dft(vel2(itheta,:),     ZERO*vel2(itheta,:),     vel2_real_tmp,     vel2_imag_tmp    )
                call dft(vel3(itheta,:),     ZERO*vel3(itheta,:),     vel3_real_tmp,     vel3_imag_tmp    )
                call dft(pressure(itheta,:), ZERO*pressure(itheta,:), pressure_real_tmp, pressure_imag_tmp)

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
        class(outlet_characteristic_quasi3d_unsteady_HB_t),  intent(inout)   :: self
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






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/12/2018
    !!
    !---------------------------------------------------------------------------------
    subroutine apply_nonreflecting_condition(self,worker,bc_comm,                   &
                                             density_Fts_real,  density_Fts_imag,   &
                                             vel1_Fts_real,     vel1_Fts_imag,      &
                                             vel2_Fts_real,     vel2_Fts_imag,      &
                                             vel3_Fts_real,     vel3_Fts_imag,      &
                                             pressure_Fts_real, pressure_Fts_imag)
        class(outlet_characteristic_quasi3d_unsteady_HB_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        type(mpi_comm),                             intent(in)      :: bc_comm
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

        integer(ik)             :: iradius, itheta, itime, ntheta
        real(rk),   allocatable :: pitch(:), spatial_periodicity(:), p_user(:)
        real(rk)                :: omega, kz
        type(AD_D)              :: state_real(5), state_imag(5), alpha_real(5), alpha_imag(5),  &
                                   density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar,     &
                                   k123, k4, k5, pyramid, c5, ddensity, dvel3, dpressure, c_bar, T(5,5), Tinv(5,5)

        p_user              = self%bcproperties%compute('Average Pressure',worker%time(),worker%coords())
        pitch               = self%bcproperties%compute('Pitch',time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])
        spatial_periodicity = self%bcproperties%compute('Spatial Periodicity', time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])

        ! Initialize derivatives for eigenvectors
        T    = ZERO*density_Fts_real(1,1,1)
        Tinv = ZERO*density_Fts_real(1,1,1)

        do iradius = 1,size(density_Fts_real,1)

            ! Get spatio-temporal average
            density_bar  = density_Fts_real(iradius,1,1)
            vel1_bar     = vel1_Fts_real(iradius,1,1)
            vel2_bar     = vel2_Fts_real(iradius,1,1)
            vel3_bar     = vel3_Fts_real(iradius,1,1)
            pressure_bar = pressure_Fts_real(iradius,1,1)
            c_bar = sqrt(gam*pressure_bar/density_bar)

            do itime = 1,size(density_Fts_real,3)

                if (itime == 1) then
                    omega = ZERO
                else
                    print*, 'WARNING: accesing frequency incorrectly'
                    !omega = worker%time_manager%freqs(itime-1)
                    omega = worker%time_manager%freqs(1)
                end if

                ntheta = size(density_Fts_real,2)
                !do itheta = 1,ntheta
                print*, 'Warning: only accounting for a single spatial mode'
                do itheta = 1,1


                    ! Handle all other modes with nonreflecting outlet condition
                    if (itime > 1) then

                        ! Compute transverse wave number 
                        ! **** WARNING: spatial_periodicity here might not be general enough! ****
                        !kz = TWO*PI*real(itheta-1,rk)/spatial_periodicity(1)
                        if (itheta <= ((ntheta-1)/2 + 1)) then
                            kz = TWO*PI*real(itheta-1,rk)/spatial_periodicity(1)  ! positive frequencies
                        else
                            kz = -TWO*PI*real(ntheta-itheta+1,rk)/spatial_periodicity(1) ! negative frequencies
                        end if
                        
                        !print*, 'omega: ', omega
                        !print*, 'kz: ', kz
                        !print*, 'c_bar: ', c_bar%x_ad_
                        !print*, 'vel3_bar: ', vel3_bar%x_ad_

                        ! Compute wave number for convected modes
                        k123 = (omega - kz*vel2_bar)/vel3_bar

                        ! Compute wave number for acoustic modes
                        ! **** WARNING: ASSUMING PYRAMID IS ALWAYS POSITIVE!!!! ****
                        pyramid = (omega - kz*vel2_bar)**TWO - kz*kz*(c_bar**TWO - vel3_bar**TWO)
                        if (pyramid < ZERO) print*, 'WARNING! pyramid < 0'
                        k4 = (-vel3_bar*(omega-kz*vel2_bar) + c_bar*sqrt(pyramid))/(c_bar**TWO - vel3_bar**TWO)
                        k5 = (-vel3_bar*(omega-kz*vel2_bar) - c_bar*sqrt(pyramid))/(c_bar**TWO - vel3_bar**TWO)
                        
                        ! Assemble right eigenvectors
                        T(1,1) = density_bar
                        T(2,1) = ZERO
                        T(3,1) = ZERO
                        T(4,1) = ZERO
                        T(5,1) = ZERO

                        T(1,2) = ZERO
                        T(2,2) = ZERO
                        T(3,2) = c_bar
                        T(4,2) = ZERO
                        T(5,2) = ZERO

                        T(1,3) = ZERO
                        T(2,3) = -c_bar*kz
                        T(3,3) = ZERO
                        T(4,3) = c_bar*k123
                        T(5,3) = ZERO

                        T(1,4) = density_bar
                        T(2,4) = -c_bar*c_bar*k4/(vel3_bar*k4 - vel3_bar*k123)
                        T(3,4) = ZERO
                        T(4,4) = -c_bar*c_bar*kz/(vel3_bar*k4 - vel3_bar*k123)
                        T(5,4) = density_bar*c_bar*c_bar

                        T(1,5) = density_bar
                        T(2,5) = -c_bar*c_bar*k5/(vel3_bar*k5 - vel3_bar*k123)
                        T(3,5) = ZERO
                        T(4,5) = -c_bar*c_bar*kz/(vel3_bar*k5 - vel3_bar*k123)
                        T(5,5) = density_bar*c_bar*c_bar


                        ! Assemble left eigenvectors
                        Tinv(1,1) = ONE/density_bar
                        Tinv(2,1) = ZERO
                        Tinv(3,1) = ZERO
                        Tinv(4,1) = ZERO
                        Tinv(5,1) = ZERO

                        Tinv(1,2) = ZERO
                        Tinv(2,2) = ZERO
                        Tinv(3,2) = -kz/(c_bar*(k123*k123 + kz*kz))
                        Tinv(4,2) = -vel3_bar*(k4*k123 - k123*k123)/(TWO*c_bar*c_bar*(k4*k123 + kz*kz))
                        Tinv(5,2) = -vel3_bar*(k5*k123 - k123*k123)/(TWO*c_bar*c_bar*(k5*k123 + kz*kz))

                        Tinv(1,3) = ZERO
                        Tinv(2,3) = ONE/c_bar
                        Tinv(3,3) = ZERO
                        Tinv(4,3) = ZERO
                        Tinv(5,3) = ZERO

                        Tinv(1,4) = ZERO
                        Tinv(2,4) = ZERO
                        Tinv(3,4) = k123/(c_bar*(k123*k123 + kz*kz))
                        Tinv(4,4) = -vel3_bar*(k4*kz - k123*kz)/(TWO*c_bar*c_bar*(k4*k123 + kz*kz))
                        Tinv(5,4) = -vel3_bar*(k5*kz - k123*kz)/(TWO*c_bar*c_bar*(k5*k123 + kz*kz))

                        Tinv(1,5) = -ONE/(density_bar*c_bar*c_bar)
                        Tinv(2,5) = ZERO
                        Tinv(3,5) = -kz/(density_bar*c_bar*vel3_bar*(k123*k123 + kz*kz))
                        Tinv(4,5) = ONE/(TWO*density_bar*c_bar*c_bar)
                        Tinv(5,5) = ONE/(TWO*density_bar*c_bar*c_bar)

                        ! Assemble state for current Fourier mode
                        state_real(1:5) = [density_Fts_real(iradius,itheta,itime),  &
                                           vel1_Fts_real(iradius,itheta,itime),     &
                                           vel2_Fts_real(iradius,itheta,itime),     &
                                           vel3_Fts_real(iradius,itheta,itime),     &
                                           pressure_Fts_real(iradius,itheta,itime)]
                        state_imag(1:5) = [density_Fts_imag(iradius,itheta,itime),  &
                                           vel1_Fts_imag(iradius,itheta,itime),     &
                                           vel2_Fts_imag(iradius,itheta,itime),     &
                                           vel3_Fts_imag(iradius,itheta,itime),     &
                                           pressure_Fts_imag(iradius,itheta,itime)]

                        ! Measure composition of interior solution in terms of the right eigenvectors by multiplying 
                        ! by the left eigenvectors
                        alpha_real = matmul(Tinv,state_real)
                        alpha_imag = matmul(Tinv,state_imag)

                        ! ***** WARNING: ASSUMING alpha(4) is always outgoing!! ******
                        ! ***** WARNING: ASSUMING alpha(5) is always incoming!! ******
                        alpha_real(5) = ZERO
                        alpha_imag(5) = ZERO

                        ! Reconstruct primitive Fourier modes from absorbing eigenvectors
                        state_real(1:5) = matmul(T,alpha_real)
                        state_imag(1:5) = matmul(T,alpha_imag)

                        !print*, 'state_real: '
                        !print*, state_real(:)%x_ad_
                        !print*, 'state_imag: '
                        !print*, state_imag(:)%x_ad_

                        ! Store
                        density_Fts_real(iradius,itheta,itime)  = state_real(1)
                        vel1_Fts_real(iradius,itheta,itime)     = state_real(2)
                        vel2_Fts_real(iradius,itheta,itime)     = state_real(3)
                        vel3_Fts_real(iradius,itheta,itime)     = state_real(4)
                        pressure_Fts_real(iradius,itheta,itime) = state_real(5)

                        density_Fts_imag(iradius,itheta,itime)  = state_imag(1)
                        vel1_Fts_imag(iradius,itheta,itime)     = state_imag(2)
                        vel2_Fts_imag(iradius,itheta,itime)     = state_imag(3)
                        vel3_Fts_imag(iradius,itheta,itime)     = state_imag(4)
                        pressure_Fts_imag(iradius,itheta,itime) = state_imag(5)

                    end if ! .not. spatio-temporal average
                end do !itheta
            end do !itime


        end do !iradius


    end subroutine apply_nonreflecting_condition
    !********************************************************************************





    !>
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_steady_nrbc(self,worker,bc_comm, &
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
        class(outlet_characteristic_quasi3d_unsteady_HB_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        type(mpi_comm),                             intent(in)      :: bc_comm
        type(AD_D),                                 intent(inout)   :: density_Fts_real(:,:)
        type(AD_D),                                 intent(inout)   :: density_Fts_imag(:,:)
        type(AD_D),                                 intent(inout)   :: vel1_Fts_real(:,:)
        type(AD_D),                                 intent(inout)   :: vel1_Fts_imag(:,:)
        type(AD_D),                                 intent(inout)   :: vel2_Fts_real(:,:)
        type(AD_D),                                 intent(inout)   :: vel2_Fts_imag(:,:)
        type(AD_D),                                 intent(inout)   :: vel3_Fts_real(:,:)
        type(AD_D),                                 intent(inout)   :: vel3_Fts_imag(:,:)
        type(AD_D),                                 intent(inout)   :: pressure_Fts_real(:,:)
        type(AD_D),                                 intent(inout)   :: pressure_Fts_imag(:,:)


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                                      &
            c1,    c2,    c3,    c4,    c5,                                             &
            c1_3d, c2_3d, c3_3d, c4_3d, c5_3d,                                          &
            c1_1d, c2_1d, c3_1d, c4_1d, c5_1d,                                          &
            c_bar, ddensity, dvel1, dvel2, dvel3, dpressure, expect_zero,               &
            density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar

        type(AD_D), allocatable, dimension(:,:) ::                              &
            c1_hat_real, c2_hat_real, c3_hat_real, c4_hat_real, c5_hat_real,    &
            c1_hat_imag, c2_hat_imag, c3_hat_imag, c4_hat_imag, c5_hat_imag

        type(AD_D)  :: pressure_avg, vel1_avg, vel2_avg, vel3_avg, density_avg, c_avg,              &
                       density_bar_r, vel1_bar_r, vel2_bar_r, vel3_bar_r, pressure_bar_r, c_bar_r,  &
                       A3_real, A3_imag, A4_real, A4_imag, beta

        type(point_t),  allocatable                 :: coords(:)
        real(rk),       allocatable, dimension(:)   :: p_user, r, pitch
        real(rk)                                    :: theta_offset
        integer(ik)                                 :: iradius, igq, ierr, imode, nmodes, itheta

        ! Get spatio-temporal average at radial stations
        density_bar  = density_Fts_real(:,1)
        vel1_bar     = vel1_Fts_real(:,1)
        vel2_bar     = vel2_Fts_real(:,1)
        vel3_bar     = vel3_Fts_real(:,1)
        pressure_bar = pressure_Fts_real(:,1)

        ! Retrieve target average pressure
        p_user = self%bcproperties%compute('Average Pressure',worker%time(),worker%coords())
        pitch  = self%bcproperties%compute('Pitch',worker%time(),worker%coords())


        ! Compute Fourier decomposition at set of radial stations: 
        !   : U_hat(nmodes,nradius)
        call self%compute_steady_decomposition(worker,bc_COMM,                          &
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

        ! Solve for c5 using nonreflecting condition
        nmodes = size(c5_hat_real,1)
        do iradius = 1,size(self%r)
            ! Get average parts
            density_bar_r  = density_bar(iradius)
            vel1_bar_r     = vel1_bar(iradius)
            vel2_bar_r     = vel2_bar(iradius)
            vel3_bar_r     = vel3_bar(iradius)
            pressure_bar_r = pressure_bar(iradius)
            c_bar_r        = sqrt(gam*pressure_bar_r/density_bar_r)

            ! starting with 2 here because the first mode is treated with 1D characteristics
            do imode = 2,nmodes 
                ! Account for sign(mode) in the calculation of beta. The second half of the
                ! modes are negative frequencies.
                if (imode <= (nmodes-1)/2 + 1) then
                    beta = sqrt(c_bar_r*c_bar_r  -  (vel3_bar_r*vel3_bar_r + vel2_bar_r*vel2_bar_r))
                else if (imode > (nmodes-1)/2 + 1) then
                    beta = -sqrt(c_bar_r*c_bar_r  -  (vel3_bar_r*vel3_bar_r + vel2_bar_r*vel2_bar_r))
                end if

                ! The imaginary part of beta has already been accounted for in
                ! the expressions for A2 and A3
                A3_real = -TWO*vel3_bar_r*vel2_bar_r/(vel2_bar_r*vel2_bar_r + beta*beta)
                A3_imag = -TWO*beta*vel3_bar_r/(vel2_bar_r*vel2_bar_r + beta*beta)

                A4_real = (beta*beta - vel2_bar_r*vel2_bar_r)/(beta*beta + vel2_bar_r*vel2_bar_r)
                A4_imag = -TWO*beta*vel2_bar_r/(beta*beta + vel2_bar_r*vel2_bar_r)

                c5_hat_real(imode,iradius) = (A3_real*c3_hat_real(imode,iradius) - A3_imag*c3_hat_imag(imode,iradius))  &   ! A3*c3 (real)
                                           - (A4_real*c4_hat_real(imode,iradius) - A4_imag*c4_hat_imag(imode,iradius))      ! A4*c4 (real)
                c5_hat_imag(imode,iradius) = (A3_imag*c3_hat_real(imode,iradius) + A3_real*c3_hat_imag(imode,iradius))  &   ! A3*c3 (imag)
                                           - (A4_imag*c4_hat_real(imode,iradius) + A4_real*c4_hat_imag(imode,iradius))      ! A4*c4 (imag)
            end do !imode
        end do !iradius



        ! Compute 1-4 characteristics from extrapolation: difference in radius-local mean and boundary average
        call self%compute_boundary_average(worker,bc_comm,density_bar,vel1_bar,vel2_bar,vel3_bar,pressure_bar, &
                                                          density_avg,vel1_avg,vel2_avg,vel3_avg,pressure_avg)
        c_avg = sqrt(gam*pressure_avg/density_avg)


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
        ! Add 1D contribution to Fourier modes
        do iradius = 1,size(self%r)
            c1_1d(iradius) = -c_avg*c_avg*ddensity(iradius)    +  dpressure(iradius)
            c2_1d(iradius) = density_avg*c_avg*dvel1(iradius)
            c3_1d(iradius) = density_avg*c_avg*dvel2(iradius)
            c4_1d(iradius) = density_avg*c_avg*dvel3(iradius)  +  dpressure(iradius)
            c5_1d(iradius) = -TWO*(pressure_avg - p_user(1))
        end do
        c1_hat_real(1,:) = c1_hat_real(1,:) + c1_1d(:)
        c2_hat_real(1,:) = c2_hat_real(1,:) + c2_1d(:)
        c3_hat_real(1,:) = c3_hat_real(1,:) + c3_1d(:)
        c4_hat_real(1,:) = c4_hat_real(1,:) + c4_1d(:)
        c5_hat_real(1,:) = c5_hat_real(1,:) + c5_1d(:)


        ! Convert characteristic Fourier modes back to primitive Fourier modes and store
        do iradius = 1,size(self%r)
            ! Get average parts
            density_bar_r  = density_bar(iradius)
            vel1_bar_r     = vel1_bar(iradius)
            vel2_bar_r     = vel2_bar(iradius)
            vel3_bar_r     = vel3_bar(iradius)
            pressure_bar_r = pressure_bar(iradius)
            c_bar_r        = sqrt(gam*pressure_bar_r/density_bar_r)
            do itheta = 1,size(c1_hat_real,1)
                density_Fts_real(iradius,itheta) = -(ONE/(c_bar_r*c_bar_r))*c1_hat_real(itheta,iradius) + (ONE/(TWO*c_bar_r*c_bar_r))*c4_hat_real(itheta,iradius) + (ONE/(TWO*c_bar_r*c_bar_r))*c5_hat_real(itheta,iradius)
                vel1_Fts_real(iradius,itheta) = (ONE/(density_bar_r*c_bar_r))*c2_hat_real(itheta,iradius)
                vel2_Fts_real(iradius,itheta) = (ONE/(density_bar_r*c_bar_r))*c3_hat_real(itheta,iradius)
                vel3_Fts_real(iradius,itheta) = (ONE/(TWO*density_bar_r*c_bar_r))*c4_hat_real(itheta,iradius) - (ONE/(TWO*density_bar_r*c_bar_r))*c5_hat_real(itheta,iradius)
                pressure_Fts_real(iradius,itheta) = HALF*c4_hat_real(itheta,iradius) + HALF*c5_hat_real(itheta,iradius)

                density_Fts_imag(iradius,itheta) = -(ONE/(c_bar_r*c_bar_r))*c1_hat_imag(itheta,iradius) + (ONE/(TWO*c_bar_r*c_bar_r))*c4_hat_imag(itheta,iradius) + (ONE/(TWO*c_bar_r*c_bar_r))*c5_hat_imag(itheta,iradius)
                vel1_Fts_imag(iradius,itheta) = (ONE/(density_bar_r*c_bar_r))*c2_hat_imag(itheta,iradius)
                vel2_Fts_imag(iradius,itheta) = (ONE/(density_bar_r*c_bar_r))*c3_hat_imag(itheta,iradius)
                vel3_Fts_imag(iradius,itheta) = (ONE/(TWO*density_bar_r*c_bar_r))*c4_hat_imag(itheta,iradius) - (ONE/(TWO*density_bar_r*c_bar_r))*c5_hat_imag(itheta,iradius)
                pressure_Fts_imag(iradius,itheta) = HALF*c4_hat_imag(itheta,iradius) + HALF*c5_hat_imag(itheta,iradius)
            end do
        end do 

    end subroutine compute_steady_nrbc
    !********************************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_steady_decomposition(self,worker,bc_comm,                    &
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
        class(outlet_characteristic_quasi3d_unsteady_HB_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        type(mpi_comm),                             intent(in)      :: bc_comm
        type(AD_D),                                 intent(in)      :: density_Fts_real(:,:)
        type(AD_D),                                 intent(in)      :: density_Fts_imag(:,:)
        type(AD_D),                                 intent(in)      :: vel1_Fts_real(:,:)
        type(AD_D),                                 intent(in)      :: vel1_Fts_imag(:,:)
        type(AD_D),                                 intent(in)      :: vel2_Fts_real(:,:)
        type(AD_D),                                 intent(in)      :: vel2_Fts_imag(:,:)
        type(AD_D),                                 intent(in)      :: vel3_Fts_real(:,:)
        type(AD_D),                                 intent(in)      :: vel3_Fts_imag(:,:)
        type(AD_D),                                 intent(in)      :: pressure_Fts_real(:,:)
        type(AD_D),                                 intent(in)      :: pressure_Fts_imag(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c1_hat_real(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c1_hat_imag(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c2_hat_real(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c2_hat_imag(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c3_hat_real(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c3_hat_imag(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c4_hat_real(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c4_hat_imag(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c5_hat_real(:,:)
        type(AD_D),     allocatable,                intent(inout)   :: c5_hat_imag(:,:)

        type(AD_D), allocatable,    dimension(:)    ::                                          &
            density, mom1, mom2, mom3, energy, vel1, vel2, vel3, pressure,                      &
            c1_real_tmp,      c2_real_tmp,   c3_real_tmp,   c4_real_tmp,   c5_real_tmp,         &
            c1_imag_tmp,      c2_imag_tmp,   c3_imag_tmp,   c4_imag_tmp,   c5_imag_tmp,         &
            density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar

        type(AD_D)  :: c_bar

        integer(ik)             :: nmodes, nradius, iradius, itheta, ntheta, imode, ierr
        real(rk)                :: shift_r, shift_i
        real(rk),   allocatable :: pitch(:)

        ! Get spatio-temporal average at radial stations
        density_bar  = density_Fts_real(:,1)
        vel1_bar     = vel1_Fts_real(:,1)
        vel2_bar     = vel2_Fts_real(:,1)
        vel3_bar     = vel3_Fts_real(:,1)
        pressure_bar = pressure_Fts_real(:,1)

        ! Define Fourier discretization
        nmodes  = self%nfourier_space
        ntheta  = 1 + (nmodes-1)*2
        nradius = size(self%r)

        pitch  = self%bcproperties%compute('Pitch',time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])

        ! Allocate storage in result
        allocate(c1_hat_real(ntheta,nradius), c1_hat_imag(ntheta,nradius),  &
                 c2_hat_real(ntheta,nradius), c2_hat_imag(ntheta,nradius),  &
                 c3_hat_real(ntheta,nradius), c3_hat_imag(ntheta,nradius),  &
                 c4_hat_real(ntheta,nradius), c4_hat_imag(ntheta,nradius),  &
                 c5_hat_real(ntheta,nradius), c5_hat_imag(ntheta,nradius), stat=ierr)
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
                c1_hat_real(itheta,iradius) = -(c_bar*c_bar)*density_Fts_real(iradius,itheta)             +  (ONE)*pressure_Fts_real(iradius,itheta)
                c2_hat_real(itheta,iradius) = (density_bar(iradius)*c_bar)*vel1_Fts_real(iradius,itheta)
                c3_hat_real(itheta,iradius) = (density_bar(iradius)*c_bar)*vel2_Fts_real(iradius,itheta)
                c4_hat_real(itheta,iradius) = (density_bar(iradius)*c_bar)*vel3_Fts_real(iradius,itheta)  +  (ONE)*pressure_Fts_real(iradius,itheta)
                c5_hat_real(itheta,iradius) = -(density_bar(iradius)*c_bar)*vel3_Fts_real(iradius,itheta) +  (ONE)*pressure_Fts_real(iradius,itheta)

                c1_hat_imag(itheta,iradius) = -(c_bar*c_bar)*density_Fts_imag(iradius,itheta)             +  (ONE)*pressure_Fts_imag(iradius,itheta)
                c2_hat_imag(itheta,iradius) = (density_bar(iradius)*c_bar)*vel1_Fts_imag(iradius,itheta)
                c3_hat_imag(itheta,iradius) = (density_bar(iradius)*c_bar)*vel2_Fts_imag(iradius,itheta)
                c4_hat_imag(itheta,iradius) = (density_bar(iradius)*c_bar)*vel3_Fts_imag(iradius,itheta)  +  (ONE)*pressure_Fts_imag(iradius,itheta)
                c5_hat_imag(itheta,iradius) = -(density_bar(iradius)*c_bar)*vel3_Fts_imag(iradius,itheta) +  (ONE)*pressure_Fts_imag(iradius,itheta)
            end do
        end do !iradius

    end subroutine compute_steady_decomposition
    !*********************************************************************************

!    !>
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   2/25/2018
!    !!
!    !--------------------------------------------------------------------------------
!    subroutine compute_steady_decomposition(self,worker,bc_comm,                                        &
!                                            density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar,    &
!                                            c1_hat_real,       c1_hat_imag,                             &
!                                            c2_hat_real,       c2_hat_imag,                             &
!                                            c3_hat_real,       c3_hat_imag,                             &
!                                            c4_hat_real,       c4_hat_imag,                             &
!                                            c5_hat_real,       c5_hat_imag)
!        class(outlet_characteristic_quasi3d_unsteady_HB_t),  intent(inout)   :: self
!        type(chidg_worker_t),                       intent(inout)   :: worker
!        type(mpi_comm),                             intent(in)      :: bc_comm
!        type(AD_D),                                 intent(in)      :: density_bar(:)
!        type(AD_D),                                 intent(in)      :: vel1_bar(:)
!        type(AD_D),                                 intent(in)      :: vel2_bar(:)
!        type(AD_D),                                 intent(in)      :: vel3_bar(:)
!        type(AD_D),                                 intent(in)      :: pressure_bar(:)
!        type(AD_D),     allocatable,                intent(inout)   :: c1_hat_real(:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: c1_hat_imag(:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: c2_hat_real(:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: c2_hat_imag(:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: c3_hat_real(:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: c3_hat_imag(:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: c4_hat_real(:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: c4_hat_imag(:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: c5_hat_real(:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: c5_hat_imag(:,:)
!
!        type(AD_D), allocatable,    dimension(:)    ::                                          &
!            density, mom1, mom2, mom3, energy, vel1, vel2, vel3, pressure,                      &
!            c1_real_tmp,      c2_real_tmp,   c3_real_tmp,   c4_real_tmp,   c5_real_tmp,         &
!            c1_imag_tmp,      c2_imag_tmp,   c3_imag_tmp,   c4_imag_tmp,   c5_imag_tmp,         &
!            c1,         c2,         c3,         c4,         c5,                                 &
!            ddensity,   dvel1,      dvel2,      dvel3,      dpressure
!
!        type(AD_D)  :: c_bar
!
!        integer(ik)             :: nmodes, nradius, iradius, itheta, ntheta, imode, ierr
!        real(rk)                :: shift_r, shift_i
!        real(rk),   allocatable :: pitch(:)
!
!        ! Define Fourier discretization
!        nmodes  = self%nfourier_space
!        ntheta  = 1 + (nmodes-1)*2
!        nradius = size(self%r)
!
!        pitch  = self%bcproperties%compute('Pitch',time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])
!
!        ! Allocate storage in result
!        allocate(c1_hat_real(ntheta,nradius), c1_hat_imag(ntheta,nradius),  &
!                 c2_hat_real(ntheta,nradius), c2_hat_imag(ntheta,nradius),  &
!                 c3_hat_real(ntheta,nradius), c3_hat_imag(ntheta,nradius),  &
!                 c4_hat_real(ntheta,nradius), c4_hat_imag(ntheta,nradius),  &
!                 c5_hat_real(ntheta,nradius), c5_hat_imag(ntheta,nradius), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!        ! Perform Fourier decomposition at each radial station.
!        do iradius = 1,nradius
!
!            ! Compute spatio-temporal average speed of sound at current radius
!            c_bar = sqrt(gam*pressure_bar(iradius)/density_bar(iradius))
!
!            ! Interpolate solution to physical_nodes at current radial station
!            density = worker%interpolate_field_general('Density',    donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:))
!            mom1    = worker%interpolate_field_general('Momentum-1', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:))
!            mom2    = worker%interpolate_field_general('Momentum-2', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:))
!            mom3    = worker%interpolate_field_general('Momentum-3', donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:))
!            energy  = worker%interpolate_field_general('Energy',     donors=self%donor(iradius,:), donor_nodes=self%donor_node(iradius,:,:))
!
!            if (worker%coordinate_system() == 'Cylindrical') then
!                mom2 = mom2/self%r(iradius)  ! convert to tangential momentum
!            end if
!
!            ! Compute velocities and pressure
!            vel1 = mom1/density
!            vel2 = mom2/density
!            vel3 = mom3/density
!            pressure = (gam-ONE)*(energy - HALF*( mom1*mom1 + mom2*mom2 + mom3*mom3 )/density )
!
!            ! Compute perturbation
!            ddensity  = density  - density_bar(iradius)
!            dvel1     = vel1     - vel1_bar(iradius)
!            dvel2     = vel2     - vel2_bar(iradius)
!            dvel3     = vel3     - vel3_bar(iradius)
!            dpressure = pressure - pressure_bar(iradius)
!
!            ! Convert perturbation to 1D characteristics
!            c1 = ZERO*ddensity
!            c2 = ZERO*ddensity
!            c3 = ZERO*ddensity
!            c4 = ZERO*ddensity
!            c5 = ZERO*ddensity
!            do itheta = 1,size(ddensity)
!                c1(itheta) = -(c_bar*c_bar)*ddensity(itheta)             +  (ONE)*dpressure(itheta)
!                c2(itheta) = (density_bar(iradius)*c_bar)*dvel1(itheta)
!                c3(itheta) = (density_bar(iradius)*c_bar)*dvel2(itheta)
!                c4(itheta) = (density_bar(iradius)*c_bar)*dvel3(itheta)  +  (ONE)*dpressure(itheta)
!                c5(itheta) = -(density_bar(iradius)*c_bar)*dvel3(itheta) +  (ONE)*dpressure(itheta)
!            end do
!
!            ! Compute Fourier transform of characteristic variables
!            call dft(c1, ZERO*c1, c1_real_tmp, c1_imag_tmp)
!            call dft(c2, ZERO*c2, c2_real_tmp, c2_imag_tmp)
!            call dft(c3, ZERO*c3, c3_real_tmp, c3_imag_tmp)
!            call dft(c4, ZERO*c4, c4_real_tmp, c4_imag_tmp)
!            call dft(c5, ZERO*c5, c5_real_tmp, c5_imag_tmp)
!
!
!            ! Adjust Fourier coefficients so their phase is relative to self%theta_ref
!            ! instead of the minimum theta of the transform.
!            !
!            !       q(relative to theta_ref) = q(relative to theta_min) * e^(j 2pi imode delta_theta/pitch)
!            !
!            ! NOTE: self%theta(:,1) are defined to be the DFT-theta_min at each radius
!            !
!            do imode = 1,size(c1_hat_real,1)
!                shift_r = realpart(exp(cmplx(ZERO,ONE)*real(imode-1,rk)*TWO*PI*(self%theta_ref-self%theta(iradius,1))/pitch(1)))
!                shift_i = imagpart(exp(cmplx(ZERO,ONE)*real(imode-1,rk)*TWO*PI*(self%theta_ref-self%theta(iradius,1))/pitch(1)))
!
!                c1_hat_real(imode,iradius) = c1_real_tmp(imode)*shift_r - c1_imag_tmp(imode)*shift_i
!                c2_hat_real(imode,iradius) = c2_real_tmp(imode)*shift_r - c2_imag_tmp(imode)*shift_i
!                c3_hat_real(imode,iradius) = c3_real_tmp(imode)*shift_r - c3_imag_tmp(imode)*shift_i
!                c4_hat_real(imode,iradius) = c4_real_tmp(imode)*shift_r - c4_imag_tmp(imode)*shift_i
!                c5_hat_real(imode,iradius) = c5_real_tmp(imode)*shift_r - c5_imag_tmp(imode)*shift_i
!
!                c1_hat_imag(imode,iradius) = c1_imag_tmp(imode)*shift_r + c1_real_tmp(imode)*shift_i
!                c2_hat_imag(imode,iradius) = c2_imag_tmp(imode)*shift_r + c2_real_tmp(imode)*shift_i
!                c3_hat_imag(imode,iradius) = c3_imag_tmp(imode)*shift_r + c3_real_tmp(imode)*shift_i
!                c4_hat_imag(imode,iradius) = c4_imag_tmp(imode)*shift_r + c4_real_tmp(imode)*shift_i
!                c5_hat_imag(imode,iradius) = c5_imag_tmp(imode)*shift_r + c5_real_tmp(imode)*shift_i
!            end do !imode
!
!        end do !iradius
!
!    end subroutine compute_steady_decomposition
!    !*********************************************************************************





    !> Compute boundary average by averaging spatio-temporal average over
    !! radial stations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/16/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_boundary_average(self,worker,bc_comm,density_bar,vel1_bar,vel2_bar,vel3_bar,pressure_bar, &
                                                            density_avg,vel1_avg,vel2_avg,vel3_avg,pressure_avg)
        class(outlet_characteristic_quasi3d_unsteady_HB_t),  intent(inout)   :: self
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
        class(outlet_characteristic_quasi3d_unsteady_HB_t),  intent(inout)   :: self
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
                    call chidg_signal(FATAL,"outlet_characteristic_quasi3d_unsteady_HB: analyze_bc_geometry, invalid face indec.")
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
        class(outlet_characteristic_quasi3d_unsteady_HB_t),  intent(inout)   :: self
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
                                   face_info_constructor(0,0,0,0,0,NO_PROC),        &   ! we don't really have a receiver face
                                   self%donor(iradius,itheta),              &
                                   self%donor_node(iradius,itheta,1:3),     &
                                   donor_found)

                ! Try LOCAL elements with try_offset if still not found 
                if ( .not. donor_found ) then
                    call find_gq_donor(mesh,                                    &
                                       node,                                    &
                                       try_offset,                              &
                                       face_info_constructor(0,0,0,0,0,NO_PROC),        &   ! we don't really have a receiver face
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
                                                face_info_constructor(0,0,0,0,0,NO_PROC),       &   ! we don't really have a receiver face
                                                self%donor(iradius,itheta),             &
                                                self%donor_node(iradius,itheta,1:3),    &
                                                donor_found)
                end if

                
                ! Try PARALLEL_ELEMENTS with try_offset if still not found 
                if ( .not. donor_found ) then
                    call find_gq_donor_parallel(mesh,                                   &
                                                node,                                   &
                                                try_offset,                             &
                                                face_info_constructor(0,0,0,0,0,NO_PROC),       &   ! we don't really have a receiver face
                                                self%donor(iradius,itheta),             &
                                                self%donor_node(iradius,itheta,1:3),    &
                                                donor_found)
                    if (donor_found) then
                        noverset=noverset+1
                        self%theta(iradius,itheta) = self%theta(iradius,itheta) - pitch(1)
                    end if
                end if 


                ! Abort if we didn't find a donor
                user_msg = "bc_state_outlet_characteristic_quasi3d_unsteady_HB%initialize_fourier_discretization: &
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




end module bc_state_outlet_characteristic_quasi3d_unsteady_HB
