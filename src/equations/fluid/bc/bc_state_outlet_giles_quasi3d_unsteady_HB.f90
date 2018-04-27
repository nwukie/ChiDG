module bc_state_outlet_giles_quasi3d_unsteady_HB
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF, ME, CYLINDRICAL,    &
                                      XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, PI
    use mod_fluid,              only: gam
    use mod_interpolation,      only: interpolate_linear, interpolate_linear_ad
    use mod_gridspace,          only: linspace
    use mod_dft,                only: dft, idft_eval
    use mod_chimera,            only: find_gq_donor, find_gq_donor_parallel

    use type_point,             only: point_t
    use type_mesh,              only: mesh_t
    use type_bc_state,          only: bc_state_t
    !use bc_state_fluid_averaging,   only: bc_fluid_averaging_t
    use bc_giles_HB_base,       only: giles_HB_base_t
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
    type, public, extends(giles_HB_base_t) :: outlet_giles_quasi3d_unsteady_HB_t

    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: compute_bc_state     ! boundary condition function implementation

    end type outlet_giles_quasi3d_unsteady_HB_t
    !*********************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. average_pressure 
    !!  @date   2/8/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(outlet_giles_quasi3d_unsteady_HB_t),   intent(inout) :: self

        ! Set name, family
        call self%set_name('Outlet - Giles Quasi3D Unsteady HB')
        call self%set_family('Outlet')

        ! Add functions
        call self%bcproperties%add('Average Pressure',    'Required')
        call self%bcproperties%add('Pitch',               'Required')
        call self%bcproperties%add('Spatial Periodicity', 'Required')

    end subroutine init
    !********************************************************************************





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
        class(outlet_giles_quasi3d_unsteady_HB_t),  intent(inout)   :: self
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

        type(AD_D), allocatable, dimension(:,:) ::                                  &
            density_t_real, vel1_t_real, vel2_t_real, vel3_t_real, pressure_t_real, &
            density_t_imag, vel1_t_imag, vel2_t_imag, vel3_t_imag, pressure_t_imag


        type(AD_D), allocatable, dimension(:,:,:) ::                                            &
            density_Ft_real, vel1_Ft_real, vel2_Ft_real, vel3_Ft_real, pressure_Ft_real,        &
            density_Ft_imag, vel1_Ft_imag, vel2_Ft_imag, vel3_Ft_imag, pressure_Ft_imag,        &
            density_Fts_real, vel1_Fts_real, vel2_Fts_real, vel3_Fts_real, pressure_Fts_real,   &
            density_Fts_imag, vel1_Fts_imag, vel2_Fts_imag, vel3_Fts_imag, pressure_Fts_imag,   &
            density_Fts_real_gq, vel1_Fts_real_gq, vel2_Fts_real_gq, vel3_Fts_real_gq, pressure_Fts_real_gq,   &
            density_Fts_imag_gq, vel1_Fts_imag_gq, vel2_Fts_imag_gq, vel3_Fts_imag_gq, pressure_Fts_imag_gq,    &
            density_grid, vel1_grid, vel2_grid, vel3_grid, pressure_grid


        real(rk),       allocatable, dimension(:)   :: p_user, r, pitch
        real(rk)                                    :: theta_offset
        type(point_t),  allocatable                 :: coords(:)
        integer                                     :: i, ngq, ivec, imode, itheta, itime, iradius, nmodes, ierr, igq


        ! Get back pressure from function.
        p_user = self%bcproperties%compute('Average Pressure',worker%time(),worker%coords())
        pitch  = self%bcproperties%compute('Pitch',           worker%time(),worker%coords())

        ! Interpolate interior solution to face quadrature nodes
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


        ! Get primitive variables at (radius,theta,time) grid.
        call self%get_q_interior(worker,bc_comm,  &
                                 density_grid,    &
                                 vel1_grid,       &
                                 vel2_grid,       &
                                 vel3_grid,       &
                                 pressure_grid)


        ! Compute Fourier decomposition of temporal data at points
        ! on the spatial transform grid.
        !   : U_Ft(nradius,ntheta,ntime)
        call self%compute_temporal_dft(worker,bc_comm,                      &
                                       density_grid,                        &
                                       vel1_grid,                           &
                                       vel2_grid,                           &
                                       vel3_grid,                           &
                                       pressure_grid,                       &
                                       density_Ft_real,  density_Ft_imag,   &
                                       vel1_Ft_real,     vel1_Ft_imag,      &
                                       vel2_Ft_real,     vel2_Ft_imag,      &
                                       vel3_Ft_real,     vel3_Ft_imag,      &
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


        call self%compute_absorbing_outlet(worker,bc_comm,       &
                                           density_Fts_real,     &
                                           density_Fts_imag,     &
                                           vel1_Fts_real,        &
                                           vel1_Fts_imag,        &
                                           vel2_Fts_real,        &
                                           vel2_Fts_imag,        &
                                           vel3_Fts_real,        &
                                           vel3_Fts_imag,        &
                                           pressure_Fts_real,    &
                                           pressure_Fts_imag)



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
            r = worker%coordinate('1','boundary')
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



!    !>
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   4/12/2018
!    !!
!    !---------------------------------------------------------------------------------
!    subroutine apply_nonreflecting_condition(self,worker,bc_comm,                   &
!                                             density_Fts_real,  density_Fts_imag,   &
!                                             vel1_Fts_real,     vel1_Fts_imag,      &
!                                             vel2_Fts_real,     vel2_Fts_imag,      &
!                                             vel3_Fts_real,     vel3_Fts_imag,      &
!                                             pressure_Fts_real, pressure_Fts_imag)
!        class(outlet_giles_quasi3d_unsteady_HB_t),  intent(inout)   :: self
!        type(chidg_worker_t),                       intent(inout)   :: worker
!        type(mpi_comm),                             intent(in)      :: bc_comm
!        type(AD_D),     allocatable,                intent(inout)   :: density_Fts_real(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: density_Fts_imag(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: vel1_Fts_real(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: vel1_Fts_imag(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: vel2_Fts_real(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: vel2_Fts_imag(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: vel3_Fts_real(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: vel3_Fts_imag(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: pressure_Fts_real(:,:,:)
!        type(AD_D),     allocatable,                intent(inout)   :: pressure_Fts_imag(:,:,:)
!
!        integer(ik)             :: iradius, itheta, itime, ntheta
!        real(rk),   allocatable :: pitch(:), spatial_periodicity(:), p_user(:)
!        real(rk)                :: omega, kz
!        type(AD_D)              :: state_real(5), state_imag(5), alpha_real(5), alpha_imag(5),  &
!                                   density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar,     &
!                                   k123, k4, k5, pyramid, c5, ddensity, dvel3, dpressure, c_bar, T(5,5), Tinv(5,5)
!
!        p_user              = self%bcproperties%compute('Average Pressure',worker%time(),worker%coords())
!        pitch               = self%bcproperties%compute('Pitch',time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])
!        spatial_periodicity = self%bcproperties%compute('Spatial Periodicity', time=ZERO,coord=[point_t(ZERO,ZERO,ZERO)])
!
!        ! Initialize derivatives for eigenvectors
!        T    = ZERO*density_Fts_real(1,1,1)
!        Tinv = ZERO*density_Fts_real(1,1,1)
!
!        do iradius = 1,size(density_Fts_real,1)
!
!            ! Get spatio-temporal average
!            density_bar  = density_Fts_real(iradius,1,1)
!            vel1_bar     = vel1_Fts_real(iradius,1,1)
!            vel2_bar     = vel2_Fts_real(iradius,1,1)
!            vel3_bar     = vel3_Fts_real(iradius,1,1)
!            pressure_bar = pressure_Fts_real(iradius,1,1)
!            c_bar = sqrt(gam*pressure_bar/density_bar)
!
!            do itime = 1,size(density_Fts_real,3)
!
!                if (itime == 1) then
!                    omega = ZERO
!                else
!                    print*, 'WARNING: accesing frequency incorrectly'
!                    !omega = worker%time_manager%freqs(itime-1)
!                    omega = worker%time_manager%freqs(1)
!                end if
!
!                ntheta = size(density_Fts_real,2)
!                !do itheta = 1,ntheta
!                print*, 'Warning: only accounting for a single spatial mode'
!                do itheta = 1,1
!
!
!                    ! Handle all other modes with nonreflecting outlet condition
!                    if (itime > 1) then
!
!                        ! Compute transverse wave number 
!                        ! **** WARNING: spatial_periodicity here might not be general enough! ****
!                        !kz = TWO*PI*real(itheta-1,rk)/spatial_periodicity(1)
!                        if (itheta <= ((ntheta-1)/2 + 1)) then
!                            kz = TWO*PI*real(itheta-1,rk)/spatial_periodicity(1)  ! positive frequencies
!                        else
!                            kz = -TWO*PI*real(ntheta-itheta+1,rk)/spatial_periodicity(1) ! negative frequencies
!                        end if
!                        
!                        !print*, 'omega: ', omega
!                        !print*, 'kz: ', kz
!                        !print*, 'c_bar: ', c_bar%x_ad_
!                        !print*, 'vel3_bar: ', vel3_bar%x_ad_
!
!                        ! Compute wave number for convected modes
!                        k123 = (omega - kz*vel2_bar)/vel3_bar
!
!                        ! Compute wave number for acoustic modes
!                        ! **** WARNING: ASSUMING PYRAMID IS ALWAYS POSITIVE!!!! ****
!                        pyramid = (omega - kz*vel2_bar)**TWO - kz*kz*(c_bar**TWO - vel3_bar**TWO)
!                        if (pyramid < ZERO) print*, 'WARNING! pyramid < 0'
!                        k4 = (-vel3_bar*(omega-kz*vel2_bar) + c_bar*sqrt(pyramid))/(c_bar**TWO - vel3_bar**TWO)
!                        k5 = (-vel3_bar*(omega-kz*vel2_bar) - c_bar*sqrt(pyramid))/(c_bar**TWO - vel3_bar**TWO)
!                        
!                        ! Assemble right eigenvectors
!                        T(1,1) = density_bar
!                        T(2,1) = ZERO
!                        T(3,1) = ZERO
!                        T(4,1) = ZERO
!                        T(5,1) = ZERO
!
!                        T(1,2) = ZERO
!                        T(2,2) = ZERO
!                        T(3,2) = c_bar
!                        T(4,2) = ZERO
!                        T(5,2) = ZERO
!
!                        T(1,3) = ZERO
!                        T(2,3) = -c_bar*kz
!                        T(3,3) = ZERO
!                        T(4,3) = c_bar*k123
!                        T(5,3) = ZERO
!
!                        T(1,4) = density_bar
!                        T(2,4) = -c_bar*c_bar*k4/(vel3_bar*k4 - vel3_bar*k123)
!                        T(3,4) = ZERO
!                        T(4,4) = -c_bar*c_bar*kz/(vel3_bar*k4 - vel3_bar*k123)
!                        T(5,4) = density_bar*c_bar*c_bar
!
!                        T(1,5) = density_bar
!                        T(2,5) = -c_bar*c_bar*k5/(vel3_bar*k5 - vel3_bar*k123)
!                        T(3,5) = ZERO
!                        T(4,5) = -c_bar*c_bar*kz/(vel3_bar*k5 - vel3_bar*k123)
!                        T(5,5) = density_bar*c_bar*c_bar
!
!
!                        ! Assemble left eigenvectors
!                        Tinv(1,1) = ONE/density_bar
!                        Tinv(2,1) = ZERO
!                        Tinv(3,1) = ZERO
!                        Tinv(4,1) = ZERO
!                        Tinv(5,1) = ZERO
!
!                        Tinv(1,2) = ZERO
!                        Tinv(2,2) = ZERO
!                        Tinv(3,2) = -kz/(c_bar*(k123*k123 + kz*kz))
!                        Tinv(4,2) = -vel3_bar*(k4*k123 - k123*k123)/(TWO*c_bar*c_bar*(k4*k123 + kz*kz))
!                        Tinv(5,2) = -vel3_bar*(k5*k123 - k123*k123)/(TWO*c_bar*c_bar*(k5*k123 + kz*kz))
!
!                        Tinv(1,3) = ZERO
!                        Tinv(2,3) = ONE/c_bar
!                        Tinv(3,3) = ZERO
!                        Tinv(4,3) = ZERO
!                        Tinv(5,3) = ZERO
!
!                        Tinv(1,4) = ZERO
!                        Tinv(2,4) = ZERO
!                        Tinv(3,4) = k123/(c_bar*(k123*k123 + kz*kz))
!                        Tinv(4,4) = -vel3_bar*(k4*kz - k123*kz)/(TWO*c_bar*c_bar*(k4*k123 + kz*kz))
!                        Tinv(5,4) = -vel3_bar*(k5*kz - k123*kz)/(TWO*c_bar*c_bar*(k5*k123 + kz*kz))
!
!                        Tinv(1,5) = -ONE/(density_bar*c_bar*c_bar)
!                        Tinv(2,5) = ZERO
!                        Tinv(3,5) = -kz/(density_bar*c_bar*vel3_bar*(k123*k123 + kz*kz))
!                        Tinv(4,5) = ONE/(TWO*density_bar*c_bar*c_bar)
!                        Tinv(5,5) = ONE/(TWO*density_bar*c_bar*c_bar)
!
!                        ! Assemble state for current Fourier mode
!                        state_real(1:5) = [density_Fts_real(iradius,itheta,itime),  &
!                                           vel1_Fts_real(iradius,itheta,itime),     &
!                                           vel2_Fts_real(iradius,itheta,itime),     &
!                                           vel3_Fts_real(iradius,itheta,itime),     &
!                                           pressure_Fts_real(iradius,itheta,itime)]
!                        state_imag(1:5) = [density_Fts_imag(iradius,itheta,itime),  &
!                                           vel1_Fts_imag(iradius,itheta,itime),     &
!                                           vel2_Fts_imag(iradius,itheta,itime),     &
!                                           vel3_Fts_imag(iradius,itheta,itime),     &
!                                           pressure_Fts_imag(iradius,itheta,itime)]
!
!                        ! Measure composition of interior solution in terms of the right eigenvectors by multiplying 
!                        ! by the left eigenvectors
!                        alpha_real = matmul(Tinv,state_real)
!                        alpha_imag = matmul(Tinv,state_imag)
!
!                        ! ***** WARNING: ASSUMING alpha(4) is always outgoing!! ******
!                        ! ***** WARNING: ASSUMING alpha(5) is always incoming!! ******
!                        alpha_real(5) = ZERO
!                        alpha_imag(5) = ZERO
!
!                        ! Reconstruct primitive Fourier modes from absorbing eigenvectors
!                        state_real(1:5) = matmul(T,alpha_real)
!                        state_imag(1:5) = matmul(T,alpha_imag)
!
!                        !print*, 'state_real: '
!                        !print*, state_real(:)%x_ad_
!                        !print*, 'state_imag: '
!                        !print*, state_imag(:)%x_ad_
!
!                        ! Store
!                        density_Fts_real(iradius,itheta,itime)  = state_real(1)
!                        vel1_Fts_real(iradius,itheta,itime)     = state_real(2)
!                        vel2_Fts_real(iradius,itheta,itime)     = state_real(3)
!                        vel3_Fts_real(iradius,itheta,itime)     = state_real(4)
!                        pressure_Fts_real(iradius,itheta,itime) = state_real(5)
!
!                        density_Fts_imag(iradius,itheta,itime)  = state_imag(1)
!                        vel1_Fts_imag(iradius,itheta,itime)     = state_imag(2)
!                        vel2_Fts_imag(iradius,itheta,itime)     = state_imag(3)
!                        vel3_Fts_imag(iradius,itheta,itime)     = state_imag(4)
!                        pressure_Fts_imag(iradius,itheta,itime) = state_imag(5)
!
!                    end if ! .not. spatio-temporal average
!                end do !itheta
!            end do !itime
!
!
!        end do !iradius
!
!
!    end subroutine apply_nonreflecting_condition
!    !********************************************************************************



end module bc_state_outlet_giles_quasi3d_unsteady_HB
