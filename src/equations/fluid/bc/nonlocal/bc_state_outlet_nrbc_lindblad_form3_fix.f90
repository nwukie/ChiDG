module bc_state_outlet_nrbc_lindblad_form3_fix
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF, ME, CYLINDRICAL,    &
                                      XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, PI
    use mod_fluid,              only: gam
    use mod_interpolation,      only: interpolate_linear, interpolate_linear_ad
    use mod_gridspace,          only: linspace
    use mod_dft,                only: idft_eval
    use mod_chimera,            only: find_gq_donor, find_gq_donor_parallel

    use type_point,             only: point_t
    use type_mesh,              only: mesh_t
    use type_bc_state,          only: bc_state_t
    use bc_nonlocal_nrbc_lindblad_base_fix,       only: nonlocal_nrbc_lindblad_base_fix_t
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
    type, public, extends(nonlocal_nrbc_lindblad_base_fix_t) :: outlet_nrbc_lindblad_form3_fix_t

    contains

        procedure   :: init                 ! Set-up bc state with options/name etc.
        procedure   :: compute_bc_state     ! boundary condition function implementation
        procedure   :: get_q_exterior
        procedure   :: compute_absorbing_outlet
        procedure   :: apply_nonreflecting_condition
        procedure   :: compute_boundary_global_average

    end type outlet_nrbc_lindblad_form3_fix_t
    !*********************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. average_pressure 
    !!  @date   2/8/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(outlet_nrbc_lindblad_form3_fix_t),   intent(inout) :: self

        ! Set name, family
        call self%set_name('Outlet - Nonlocal NRBC Lindblad Form3 Fix')
        call self%set_family('Outlet')

        ! Add functions
        call self%bcproperties%add('Average Pressure',    'Required')
        call self%bcproperties%add('Pitch A',             'Required')
        call self%bcproperties%add('Pitch B',             'Required')
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
        class(outlet_nrbc_lindblad_form3_fix_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop
        type(mpi_comm),                             intent(in)      :: bc_comm


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                                                      &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc, pressure_bc, vel1_bc, vel2_bc, vel3_bc,   &
            density_m, mom1_m, mom2_m, mom3_m, energy_m, pressure_m, vel1_m, vel2_m, vel3_m,   &
            density_reflection, vel1_reflection, vel2_reflection, vel3_reflection, pressure_reflection, &
            c1_1d, c2_1d, c3_1d, c4_1d, c5_1d, density_1d_out, vel1_1d_out, vel2_1d_out, vel3_1d_out, pressure_1d_out, &
            ddensity, dvel1, dvel2, dvel3, dpressure, &
            density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, c_bar, &
            ddensity_1d, dvel1_1d, dvel2_1d, dvel3_1d, dpressure_1d, &
            ddensity_3d, dvel1_3d, dvel2_3d, dvel3_3d, dpressure_3d, &
            grad1_density_m, grad1_mom1_m, grad1_mom2_m, grad1_mom3_m, grad1_energy_m,                  &
            grad2_density_m, grad2_mom1_m, grad2_mom2_m, grad2_mom3_m, grad2_energy_m,                  &
            grad3_density_m, grad3_mom1_m, grad3_mom2_m, grad3_mom3_m, grad3_energy_m,                  &
            density_unresolved, vel1_unresolved, vel2_unresolved, vel3_unresolved, pressure_unresolved, &
            density_fourier_gq, vel1_fourier_gq, vel2_fourier_gq, vel3_fourier_gq, pressure_fourier_gq


        type(AD_D), allocatable, dimension(:,:) ::                                  &
            density_check_real_gq, vel1_check_real_gq, vel2_check_real_gq, vel3_check_real_gq, pressure_check_real_gq, &
            density_check_imag_gq, vel1_check_imag_gq, vel2_check_imag_gq, vel3_check_imag_gq, pressure_check_imag_gq

        type(AD_D), allocatable, dimension(:,:,:) ::                                                                                &
            density_check_real_m, vel1_check_real_m, vel2_check_real_m, vel3_check_real_m, pressure_check_real_m, c_check_real_m,   &
            density_check_imag_m, vel1_check_imag_m, vel2_check_imag_m, vel3_check_imag_m, pressure_check_imag_m, c_check_imag_m,   &
            density_hat_real_m, vel1_hat_real_m, vel2_hat_real_m, vel3_hat_real_m, pressure_hat_real_m, c_hat_real_m,               &
            density_hat_imag_m, vel1_hat_imag_m, vel2_hat_imag_m, vel3_hat_imag_m, pressure_hat_imag_m, c_hat_imag_m,               &
            density_check_real_p, vel1_check_real_p, vel2_check_real_p, vel3_check_real_p, pressure_check_real_p, c_check_real_p,   &
            density_check_imag_p, vel1_check_imag_p, vel2_check_imag_p, vel3_check_imag_p, pressure_check_imag_p, c_check_imag_p,   &
            density_hat_real_p, vel1_hat_real_p, vel2_hat_real_p, vel3_hat_real_p, pressure_hat_real_p, c_hat_real_p,               &
            density_hat_imag_p, vel1_hat_imag_p, vel2_hat_imag_p, vel3_hat_imag_p, pressure_hat_imag_p, c_hat_imag_p,               &
            density_hat_real_abs, vel1_hat_real_abs, vel2_hat_real_abs, vel3_hat_real_abs, pressure_hat_real_abs,                   &
            density_hat_imag_abs, vel1_hat_imag_abs, vel2_hat_imag_abs, vel3_hat_imag_abs, pressure_hat_imag_abs,                   &
            density_hat_real_gq, vel1_hat_real_gq, vel2_hat_real_gq, vel3_hat_real_gq, pressure_hat_real_gq,    &
            density_hat_imag_gq, vel1_hat_imag_gq, vel2_hat_imag_gq, vel3_hat_imag_gq, pressure_hat_imag_gq,    &
            density_grid_m, vel1_grid_m, vel2_grid_m, vel3_grid_m, pressure_grid_m, c_grid_m,                   &
            density_grid_p, vel1_grid_p, vel2_grid_p, vel3_grid_p, pressure_grid_p, c_grid_p

        type(AD_D)  :: density_avg, vel1_avg, vel2_avg, vel3_avg, pressure_avg, c_avg

        real(rk),       allocatable, dimension(:)   :: p_user, r
        integer(ik) :: igq



        density_m = worker%get_field('Density',    'value', 'face interior')
        mom1_m    = worker%get_field('Momentum-1', 'value', 'face interior')
        mom2_m    = worker%get_field('Momentum-2', 'value', 'face interior')
        mom3_m    = worker%get_field('Momentum-3', 'value', 'face interior')
        energy_m  = worker%get_field('Energy',     'value', 'face interior')

        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','boundary') 
            mom2_m = mom2_m / r
        end if

        ! Compute velocities
        vel1_m = mom1_m/density_m
        vel2_m = mom2_m/density_m
        vel3_m = mom3_m/density_m

        ! Get interior pressure
        pressure_m = worker%get_field('Pressure', 'value', 'face interior')


















        ! Get back pressure from function.
        p_user = self%bcproperties%compute('Average Pressure',worker%time(),worker%coords())

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
        call self%get_q_side(worker,bc_comm,'A', &
                             density_grid_m,     &
                             vel1_grid_m,        &
                             vel2_grid_m,        &
                             vel3_grid_m,        &
                             pressure_grid_m)
        c_grid_m = sqrt(gam*pressure_grid_m/density_grid_m)


        ! Compute Fourier decomposition of temporal data at points
        ! on the spatial transform grid.
        !   q_check(r,theta,omega) = DFT(q)[time]
        call self%compute_temporal_dft(worker,bc_comm,                                  &
                                       density_grid_m,                                  &
                                       vel1_grid_m,                                     &
                                       vel2_grid_m,                                     &
                                       vel3_grid_m,                                     &
                                       pressure_grid_m,                                 &
                                       c_grid_m,                                        &
                                       density_check_real_m,  density_check_imag_m,     &
                                       vel1_check_real_m,     vel1_check_imag_m,        &
                                       vel2_check_real_m,     vel2_check_imag_m,        &
                                       vel3_check_real_m,     vel3_check_imag_m,        &
                                       pressure_check_real_m, pressure_check_imag_m,    &
                                       c_check_real_m,        c_check_imag_m)


        ! Compute Fourier decomposition in theta at set of radial 
        ! stations for each temporal mode:
        !   q_hat(r,m,omega) = DFT(q_check)[theta]
        call self%compute_spatial_dft(worker,bc_comm,'A',                           &
                                      density_check_real_m,  density_check_imag_m,  &
                                      vel1_check_real_m,     vel1_check_imag_m,     &
                                      vel2_check_real_m,     vel2_check_imag_m,     &
                                      vel3_check_real_m,     vel3_check_imag_m,     &
                                      pressure_check_real_m, pressure_check_imag_m, &
                                      c_check_real_m,        c_check_imag_m,        &
                                      density_hat_real_m,    density_hat_imag_m,    &
                                      vel1_hat_real_m,       vel1_hat_imag_m,       &
                                      vel2_hat_real_m,       vel2_hat_imag_m,       &
                                      vel3_hat_real_m,       vel3_hat_imag_m,       &
                                      pressure_hat_real_m,   pressure_hat_imag_m,   &
                                      c_hat_real_m,          c_hat_imag_m)



        !!!!! Evaluate fourier solution

        ! q_abs(r_gq) = I(q_abs(r_aux))
        call self%interpolate_raux_to_rgq(worker,bc_comm,                               &
                                          density_hat_real_m,  density_hat_imag_m,  &
                                          vel1_hat_real_m,     vel1_hat_imag_m,     &
                                          vel2_hat_real_m,     vel2_hat_imag_m,     &
                                          vel3_hat_real_m,     vel3_hat_imag_m,     &
                                          pressure_hat_real_m, pressure_hat_imag_m, &
                                          density_hat_real_gq,   density_hat_imag_gq,   &
                                          vel1_hat_real_gq,      vel1_hat_imag_gq,      &
                                          vel2_hat_real_gq,      vel2_hat_imag_gq,      &
                                          vel3_hat_real_gq,      vel3_hat_imag_gq,      &
                                          pressure_hat_real_gq,  pressure_hat_imag_gq)
        

        ! Reconstruct primitive variables at quadrature nodes from absorbing Fourier modes
        ! via inverse transform.
        !   q_check(rgq,theta,omega) = IDFT(q_hat)[m]
        call self%compute_spatial_idft_gq(worker,bc_comm,'A',                               &
                                          density_hat_real_gq,      density_hat_imag_gq,    & 
                                          vel1_hat_real_gq,         vel1_hat_imag_gq,       &
                                          vel2_hat_real_gq,         vel2_hat_imag_gq,       &
                                          vel3_hat_real_gq,         vel3_hat_imag_gq,       &
                                          pressure_hat_real_gq,     pressure_hat_imag_gq,   &
                                          density_check_real_gq,    density_check_imag_gq,  &
                                          vel1_check_real_gq,       vel1_check_imag_gq,     &
                                          vel2_check_real_gq,       vel2_check_imag_gq,     &
                                          vel3_check_real_gq,       vel3_check_imag_gq,     &
                                          pressure_check_real_gq,   pressure_check_imag_gq)


        ! q(rgq,theta,t) = IDFT(q_check)[omega]
        call self%compute_temporal_idft_gq(worker,bc_comm,                                                      &
                                           density_check_real_gq,    density_check_imag_gq,     & 
                                           vel1_check_real_gq,       vel1_check_imag_gq,        &
                                           vel2_check_real_gq,       vel2_check_imag_gq,        &
                                           vel3_check_real_gq,       vel3_check_imag_gq,        &
                                           pressure_check_real_gq,   pressure_check_imag_gq,    &
                                           density_fourier_gq, vel1_fourier_gq, vel2_fourier_gq, vel3_fourier_gq, pressure_fourier_gq)




        ! Get spatio-temporal average at radial stations
        density_bar  = density_hat_real_gq(:,1,1)
        vel1_bar     = vel1_hat_real_gq(:,1,1)
        vel2_bar     = vel2_hat_real_gq(:,1,1)
        vel3_bar     = vel3_hat_real_gq(:,1,1)
        pressure_bar = pressure_hat_real_gq(:,1,1)
        c_bar = sqrt(gam*pressure_bar/density_bar)
        !c_bar        = c_hat_real_gq(:,1,1)


!        ! Get spatio-temporal average at radial stations
!        density_bar  = density_hat_real_m(:,1,1)
!        vel1_bar     = vel1_hat_real_m(:,1,1)
!        vel2_bar     = vel2_hat_real_m(:,1,1)
!        vel3_bar     = vel3_hat_real_m(:,1,1)
!        pressure_bar = pressure_hat_real_m(:,1,1)
!        c_bar        = c_hat_real_m(:,1,1)
!
!
!        ! Compute spatio-temporal average over entire surface
!        call self%compute_boundary_average(worker,bc_comm,density_bar,vel1_bar,vel2_bar,vel3_bar,pressure_bar,c_bar, &
!                                                          density_avg,vel1_avg,vel2_avg,vel3_avg,pressure_avg,c_avg)




        density_unresolved  = density_m  - density_fourier_gq
        vel1_unresolved     = vel1_m     - vel1_fourier_gq
        vel2_unresolved     = vel2_m     - vel2_fourier_gq
        vel3_unresolved     = vel3_m     - vel3_fourier_gq
        pressure_unresolved = pressure_m - pressure_fourier_gq


        ! compute 1d characteristics associated with unresolved perturbations
        c1_1d = ZERO*density_m
        c2_1d = ZERO*density_m
        c3_1d = ZERO*density_m
        c4_1d = ZERO*density_m
        c5_1d = ZERO*density_m
!        do igq = 1,size(density_m)
!            c1_1d(igq) = -c_avg*c_avg*density_unresolved(igq)    +  pressure_unresolved(igq)
!            c2_1d(igq) = density_avg*c_avg*vel1_unresolved(igq)
!            c3_1d(igq) = density_avg*c_avg*vel2_unresolved(igq)
!            c4_1d(igq) = density_avg*c_avg*vel3_unresolved(igq)  +  pressure_unresolved(igq)
!            c5_1d(igq) = -density_avg*c_avg*vel3_unresolved(igq)  +  pressure_unresolved(igq)
!        end do
        c1_1d = -c_bar*c_bar*density_unresolved    +  pressure_unresolved
        c2_1d =  density_bar*c_bar*vel1_unresolved
        c3_1d =  density_bar*c_bar*vel2_unresolved
        c4_1d =  density_bar*c_bar*vel3_unresolved  +  pressure_unresolved
        c5_1d = -density_bar*c_bar*vel3_unresolved  +  pressure_unresolved


        ! compute 1d characteristic outgoing: c1-c4
        c5_1d = ZERO

        ! Now add local perturbation from the average
        density_1d_out  = ZERO*density_m
        vel1_1d_out     = ZERO*density_m
        vel2_1d_out     = ZERO*density_m
        vel3_1d_out     = ZERO*density_m
        pressure_1d_out = ZERO*density_m
!        do igq = 1,size(density_m)
!            density_1d_out(igq)  = (-ONE/(c_avg*c_avg))*c1_1d(igq)  +  (ONE/(TWO*c_avg*c_avg))*c4_1d(igq)  +  (ONE/(TWO*c_avg*c_avg))*c5_1d(igq)
!            vel1_1d_out(igq)     = (ONE/(density_avg*c_avg))*c2_1d(igq)
!            vel2_1d_out(igq)     = (ONE/(density_avg*c_avg))*c3_1d(igq)
!            vel3_1d_out(igq)     = (ONE/(TWO*density_avg*c_avg))*c4_1d(igq)  -  (ONE/(TWO*density_avg*c_avg))*c5_1d(igq)
!            pressure_1d_out(igq) = HALF*c4_1d(igq)  +  HALF*c5_1d(igq)
!        end do
        density_1d_out  = (-ONE/(c_bar*c_bar))*c1_1d  +  (ONE/(TWO*c_bar*c_bar))*c4_1d  +  (ONE/(TWO*c_bar*c_bar))*c5_1d
        vel1_1d_out     = (ONE/(density_bar*c_bar))*c2_1d
        vel2_1d_out     = (ONE/(density_bar*c_bar))*c3_1d
        vel3_1d_out     = (ONE/(TWO*density_bar*c_bar))*c4_1d  -  (ONE/(TWO*density_bar*c_bar))*c5_1d
        pressure_1d_out = HALF*c4_1d  +  HALF*c5_1d



















        ! Get primitive variables at (radius,theta,time) grid.
        call self%get_q_exterior(worker,bc_comm,    &
                                 density_grid_p,    &
                                 vel1_grid_p,       &
                                 vel2_grid_p,       &
                                 vel3_grid_p,       &
                                 pressure_grid_p)
        c_grid_p = sqrt(gam*pressure_grid_p/density_grid_p)


        ! Compute Fourier decomposition of temporal data at points
        ! on the spatial transform grid.
        !   q_check(r,theta,omega) = DFT(q)[time]
        call self%compute_temporal_dft(worker,bc_comm,                                  &
                                       density_grid_p,                                  &
                                       vel1_grid_p,                                     &
                                       vel2_grid_p,                                     &
                                       vel3_grid_p,                                     &
                                       pressure_grid_p,                                 &
                                       c_grid_p,                                        &
                                       density_check_real_p,  density_check_imag_p,     &
                                       vel1_check_real_p,     vel1_check_imag_p,        &
                                       vel2_check_real_p,     vel2_check_imag_p,        &
                                       vel3_check_real_p,     vel3_check_imag_p,        &
                                       pressure_check_real_p, pressure_check_imag_p,    &
                                       c_check_real_p,        c_check_imag_p)

    
        ! Compute Fourier decomposition in theta at set of radial 
        ! stations for each temporal mode:
        !   q_hat(r,m,omega) = DFT(q_check)[theta]
        call self%compute_spatial_dft(worker,bc_comm,'B',                           &
                                      density_check_real_p,  density_check_imag_p,  &
                                      vel1_check_real_p,     vel1_check_imag_p,     &
                                      vel2_check_real_p,     vel2_check_imag_p,     &
                                      vel3_check_real_p,     vel3_check_imag_p,     &
                                      pressure_check_real_p, pressure_check_imag_p, &
                                      c_check_real_p,        c_check_imag_p,        &
                                      density_hat_real_p,    density_hat_imag_p,    &
                                      vel1_hat_real_p,       vel1_hat_imag_p,       &
                                      vel2_hat_real_p,       vel2_hat_imag_p,       &
                                      vel3_hat_real_p,       vel3_hat_imag_p,       &
                                      pressure_hat_real_p,   pressure_hat_imag_p,   &
                                      c_hat_real_p,          c_hat_imag_p)


        ! Compute absorbing state using 2D eigenmode decomposition.
        !
        !   q_abs = f(q_hat_m,q_hat_p)
        !
        call self%compute_absorbing_outlet(worker,bc_comm,                              &
                                           density_hat_real_m,    density_hat_imag_m,   &
                                           vel1_hat_real_m,       vel1_hat_imag_m,      &
                                           vel2_hat_real_m,       vel2_hat_imag_m,      &
                                           vel3_hat_real_m,       vel3_hat_imag_m,      &
                                           pressure_hat_real_m,   pressure_hat_imag_m,  &
                                           c_hat_real_m,          c_hat_imag_m,         &
                                           density_hat_real_p,    density_hat_imag_p,   &
                                           vel1_hat_real_p,       vel1_hat_imag_p,      &
                                           vel2_hat_real_p,       vel2_hat_imag_p,      &
                                           vel3_hat_real_p,       vel3_hat_imag_p,      &
                                           pressure_hat_real_p,   pressure_hat_imag_p,  &
                                           c_hat_real_p,          c_hat_imag_p,         &
                                           density_hat_real_abs,  density_hat_imag_abs, &
                                           vel1_hat_real_abs,     vel1_hat_imag_abs,    &
                                           vel2_hat_real_abs,     vel2_hat_imag_abs,    &
                                           vel3_hat_real_abs,     vel3_hat_imag_abs,    &
                                           pressure_hat_real_abs, pressure_hat_imag_abs)


        ! q_abs(r_gq) = I(q_abs(r_aux))
        call self%interpolate_raux_to_rgq(worker,bc_comm,                               &
                                          density_hat_real_abs,  density_hat_imag_abs,  &
                                          vel1_hat_real_abs,     vel1_hat_imag_abs,     &
                                          vel2_hat_real_abs,     vel2_hat_imag_abs,     &
                                          vel3_hat_real_abs,     vel3_hat_imag_abs,     &
                                          pressure_hat_real_abs, pressure_hat_imag_abs, &
                                          density_hat_real_gq,   density_hat_imag_gq,   &
                                          vel1_hat_real_gq,      vel1_hat_imag_gq,      &
                                          vel2_hat_real_gq,      vel2_hat_imag_gq,      &
                                          vel3_hat_real_gq,      vel3_hat_imag_gq,      &
                                          pressure_hat_real_gq,  pressure_hat_imag_gq)


        ! Reconstruct primitive variables at quadrature nodes from absorbing Fourier modes
        ! via inverse transform.
        !   q_check(rgq,theta,omega) = IDFT(q_hat)[m]
        call self%compute_spatial_idft_gq(worker,bc_comm,'A',                               &
                                          density_hat_real_gq,      density_hat_imag_gq,    & 
                                          vel1_hat_real_gq,         vel1_hat_imag_gq,       &
                                          vel2_hat_real_gq,         vel2_hat_imag_gq,       &
                                          vel3_hat_real_gq,         vel3_hat_imag_gq,       &
                                          pressure_hat_real_gq,     pressure_hat_imag_gq,   &
                                          density_check_real_gq,    density_check_imag_gq,  &
                                          vel1_check_real_gq,       vel1_check_imag_gq,     &
                                          vel2_check_real_gq,       vel2_check_imag_gq,     &
                                          vel3_check_real_gq,       vel3_check_imag_gq,     &
                                          pressure_check_real_gq,   pressure_check_imag_gq)


        ! q(rgq,theta,t) = IDFT(q_check)[omega]
        call self%compute_temporal_idft_gq(worker,bc_comm,                                      &
                                           density_check_real_gq,    density_check_imag_gq,     & 
                                           vel1_check_real_gq,       vel1_check_imag_gq,        &
                                           vel2_check_real_gq,       vel2_check_imag_gq,        &
                                           vel3_check_real_gq,       vel3_check_imag_gq,        &
                                           pressure_check_real_gq,   pressure_check_imag_gq,    &
                                           density_bc, vel1_bc, vel2_bc, vel3_bc, pressure_bc)



        ! Add in correction from part unresolved by Fourier decomposition
        density_bc = density_bc   + density_1d_out
        vel1_bc    = vel1_bc      + vel1_1d_out
        vel2_bc    = vel2_bc      + vel2_1d_out
        vel3_bc    = vel3_bc      + vel3_1d_out
        pressure_bc = pressure_bc + pressure_1d_out





        !
        ! Form conserved variables
        !
        density_bc = density_bc
        mom1_bc    = density_bc*vel1_bc
        mom2_bc    = density_bc*vel2_bc
        mom3_bc    = density_bc*vel3_bc
        energy_bc  = pressure_bc/(gam - ONE)  + HALF*(mom1_bc*mom1_bc + mom2_bc*mom2_bc + mom3_bc*mom3_bc)/density_bc


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





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/18/2018
    !!
    !------------------------------------------------------------------------------------
    subroutine get_q_exterior(self,worker,bc_comm, density, vel1, vel2, vel3, pressure)
        class(outlet_nrbc_lindblad_form3_fix_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        type(mpi_comm),                             intent(in)      :: bc_comm
        type(AD_D),     allocatable,                intent(inout)   :: density(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel1(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel2(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: vel3(:,:,:)
        type(AD_D),     allocatable,                intent(inout)   :: pressure(:,:,:)

        type(AD_D), allocatable, dimension(:,:,:) ::  &
            mom1, mom2, mom3, energy

        real(rk),   allocatable, dimension(:)   :: p_user

        integer(ik) :: iradius, itime, nradius, ntheta, ntime, ierr


        ! Define Fourier space discretization to determine
        ! number of theta-samples are being taken
        nradius = size(self%r)
        ntheta  = size(self%theta_a,2)
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
                density(iradius,:,itime) = worker%interpolate_field_general('Density',    donors=self%donor_a(iradius,:), donor_nodes=self%donor_node_a(iradius,:,:), itime=itime)
                mom1(iradius,:,itime)    = worker%interpolate_field_general('Momentum-1', donors=self%donor_a(iradius,:), donor_nodes=self%donor_node_a(iradius,:,:), itime=itime)
                mom2(iradius,:,itime)    = worker%interpolate_field_general('Momentum-2', donors=self%donor_a(iradius,:), donor_nodes=self%donor_node_a(iradius,:,:), itime=itime)
                mom3(iradius,:,itime)    = worker%interpolate_field_general('Momentum-3', donors=self%donor_a(iradius,:), donor_nodes=self%donor_node_a(iradius,:,:), itime=itime)
                energy(iradius,:,itime)  = worker%interpolate_field_general('Energy',     donors=self%donor_a(iradius,:), donor_nodes=self%donor_node_a(iradius,:,:), itime=itime)

                if (worker%coordinate_system() == 'Cylindrical') then
                    mom2(iradius,:,itime) = mom2(iradius,:,itime)/self%r(iradius)  ! convert to tangential momentum
                end if
            end do
        end do

        ! Compute velocities and pressure at each time
        density  = density
        vel1     = mom1/density
        vel2     = mom2/density
        vel3     = mom3/density

        ! Pressure from user
        pressure = (gam-ONE)*(energy - HALF*(mom1*mom1 + mom2*mom2 + mom3*mom3)/density)
        p_user = self%bcproperties%compute('Average Pressure',worker%time(),[point_t(ZERO,ZERO,ZERO)])
        pressure = p_user(1)

    end subroutine get_q_exterior
    !************************************************************************************




    !> DANIEL'S FORMULATION
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_absorbing_outlet(self,worker,bc_comm,                    &
                                     density_real_m,    density_imag_m,      &
                                     vel1_real_m,       vel1_imag_m,         &
                                     vel2_real_m,       vel2_imag_m,         &
                                     vel3_real_m,       vel3_imag_m,         &
                                     pressure_real_m,   pressure_imag_m,     &
                                     c_real_m,          c_imag_m,            &
                                     density_real_p,    density_imag_p,      &
                                     vel1_real_p,       vel1_imag_p,         &
                                     vel2_real_p,       vel2_imag_p,         &
                                     vel3_real_p,       vel3_imag_p,         &
                                     pressure_real_p,   pressure_imag_p,     &
                                     c_real_p,          c_imag_p,            &
                                     density_real_abs,  density_imag_abs,    &
                                     vel1_real_abs,     vel1_imag_abs,       &
                                     vel2_real_abs,     vel2_imag_abs,       &
                                     vel3_real_abs,     vel3_imag_abs,       &
                                     pressure_real_abs, pressure_imag_abs)
        class(outlet_nrbc_lindblad_form3_fix_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(mpi_comm),             intent(in)      :: bc_comm
        type(AD_D),                 intent(inout)   :: density_real_m(:,:,:)
        type(AD_D),                 intent(inout)   :: density_imag_m(:,:,:)
        type(AD_D),                 intent(inout)   :: vel1_real_m(:,:,:)
        type(AD_D),                 intent(inout)   :: vel1_imag_m(:,:,:)
        type(AD_D),                 intent(inout)   :: vel2_real_m(:,:,:)
        type(AD_D),                 intent(inout)   :: vel2_imag_m(:,:,:)
        type(AD_D),                 intent(inout)   :: vel3_real_m(:,:,:)
        type(AD_D),                 intent(inout)   :: vel3_imag_m(:,:,:)
        type(AD_D),                 intent(inout)   :: pressure_real_m(:,:,:)
        type(AD_D),                 intent(inout)   :: pressure_imag_m(:,:,:)
        type(AD_D),                 intent(inout)   :: c_real_m(:,:,:)
        type(AD_D),                 intent(inout)   :: c_imag_m(:,:,:)
        type(AD_D),                 intent(inout)   :: density_real_p(:,:,:)
        type(AD_D),                 intent(inout)   :: density_imag_p(:,:,:)
        type(AD_D),                 intent(inout)   :: vel1_real_p(:,:,:)
        type(AD_D),                 intent(inout)   :: vel1_imag_p(:,:,:)
        type(AD_D),                 intent(inout)   :: vel2_real_p(:,:,:)
        type(AD_D),                 intent(inout)   :: vel2_imag_p(:,:,:)
        type(AD_D),                 intent(inout)   :: vel3_real_p(:,:,:)
        type(AD_D),                 intent(inout)   :: vel3_imag_p(:,:,:)
        type(AD_D),                 intent(inout)   :: pressure_real_p(:,:,:)
        type(AD_D),                 intent(inout)   :: pressure_imag_p(:,:,:)
        type(AD_D),                 intent(inout)   :: c_real_p(:,:,:)
        type(AD_D),                 intent(inout)   :: c_imag_p(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: density_real_abs(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: density_imag_abs(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_real_abs(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel1_imag_abs(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_real_abs(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel2_imag_abs(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_real_abs(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: vel3_imag_abs(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_real_abs(:,:,:)
        type(AD_D), allocatable,    intent(inout)   :: pressure_imag_abs(:,:,:)

        type(AD_D), allocatable, dimension(:,:,:) ::                &
            a1_real, a2_real, a3_real, a4_real, a5_real,            &
            a1_imag, a2_imag, a3_imag, a4_imag, a5_imag,            &
            a1_real_m, a2_real_m, a3_real_m, a4_real_m, a5_real_m,  &
            a1_imag_m, a2_imag_m, a3_imag_m, a4_imag_m, a5_imag_m,  &
            a1_real_p, a2_real_p, a3_real_p, a4_real_p, a5_real_p,  &
            a1_imag_p, a2_imag_p, a3_imag_p, a4_imag_p, a5_imag_p,  &
            c1_real_m, c2_real_m, c3_real_m, c4_real_m, c5_real_m,  &
            c1_imag_m, c2_imag_m, c3_imag_m, c4_imag_m, c5_imag_m


        type(AD_D)  :: density_avg, vel1_avg, vel2_avg, vel3_avg, pressure_avg, c_avg, &
                       density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, c_bar, &
                       B3_real, B3_imag, B4_real, B4_imag, beta

        type(AD_D), allocatable, dimension(:)   :: ddensity, dvel1, dvel2, dvel3, dpressure, c1_1d, c2_1d, c3_1d, c4_1d, c5_1d

        real(rk),   allocatable, dimension(:)   :: p_user
        integer(ik) :: iradius, itheta, itime, ntheta

        ! Retrieve target average pressure
        p_user = self%bcproperties%compute('Average Pressure',worker%time(),worker%coords())


        !
        ! Step 1: handle higher-order modes: m /= 0
        !
        ! Project interior to eigenmodes
        call self%primitive_to_eigenmodes(worker,bc_comm,'A',               &
                                          density_real_m(:,1,1),            &
                                          vel1_real_m(:,1,1),               &
                                          vel2_real_m(:,1,1),               &
                                          vel3_real_m(:,1,1),               &
                                          pressure_real_m(:,1,1),           &
                                          c_real_m(:,1,1),                  &
                                          density_real_m,  density_imag_m,  &
                                          vel1_real_m,     vel1_imag_m,     &
                                          vel2_real_m,     vel2_imag_m,     &
                                          vel3_real_m,     vel3_imag_m,     &
                                          pressure_real_m, pressure_imag_m, &
                                          a1_real_m,       a1_imag_m,       &
                                          a2_real_m,       a2_imag_m,       &
                                          a3_real_m,       a3_imag_m,       &
                                          a4_real_m,       a4_imag_m,       &
                                          a5_real_m,       a5_imag_m)


        ! Project exterior to eigenmodes
        call self%primitive_to_eigenmodes(worker,bc_comm,'B',               &
                                          density_real_m(:,1,1),            &
                                          vel1_real_m(:,1,1),               &
                                          vel2_real_m(:,1,1),               &
                                          vel3_real_m(:,1,1),               &
                                          pressure_real_m(:,1,1),           &
                                          c_real_m(:,1,1),                  &
                                          density_real_p,  density_imag_p,  &
                                          vel1_real_p,     vel1_imag_p,     &
                                          vel2_real_p,     vel2_imag_p,     &
                                          vel3_real_p,     vel3_imag_p,     &
                                          pressure_real_p, pressure_imag_p, &
                                          a1_real_p,       a1_imag_p,       &
                                          a2_real_p,       a2_imag_p,       &
                                          a3_real_p,       a3_imag_p,       &
                                          a4_real_p,       a4_imag_p,       &
                                          a5_real_p,       a5_imag_p)



        call self%apply_nonreflecting_condition(worker,bc_comm,             &
                                                density_real_m(:,1,1),      &
                                                vel1_real_m(:,1,1),         &
                                                vel2_real_m(:,1,1),         &
                                                vel3_real_m(:,1,1),         &
                                                pressure_real_m(:,1,1),     &
                                                c_real_m(:,1,1),            &
                                                a1_real_m,    a1_imag_m,    &
                                                a2_real_m,    a2_imag_m,    &
                                                a3_real_m,    a3_imag_m,    &
                                                a4_real_m,    a4_imag_m,    &
                                                a5_real_m,    a5_imag_m,    &
                                                a1_real_p,    a1_imag_p,    &
                                                a2_real_p,    a2_imag_p,    &
                                                a3_real_p,    a3_imag_p,    &
                                                a4_real_p,    a4_imag_p,    &
                                                a5_real_p,    a5_imag_p,    &
                                                a1_real,      a1_imag,      &
                                                a2_real,      a2_imag,      &
                                                a3_real,      a3_imag,      &
                                                a4_real,      a4_imag,      &
                                                a5_real,      a5_imag)


        ! Convert back to primitive variables
        call self%eigenmodes_to_primitive(worker,bc_comm,self%get_face_side(worker),&
                                          density_real_m(:,1,1),                &
                                          vel1_real_m(:,1,1),                   &
                                          vel2_real_m(:,1,1),                   &
                                          vel3_real_m(:,1,1),                   &
                                          pressure_real_m(:,1,1),               &
                                          c_real_m(:,1,1),                      &
                                          a1_real,           a1_imag,           &
                                          a2_real,           a2_imag,           &
                                          a3_real,           a3_imag,           &
                                          a4_real,           a4_imag,           &
                                          a5_real,           a5_imag,           &
                                          density_real_abs,  density_imag_abs,  &
                                          vel1_real_abs,     vel1_imag_abs,     &
                                          vel2_real_abs,     vel2_imag_abs,     &
                                          vel3_real_abs,     vel3_imag_abs,     &
                                          pressure_real_abs, pressure_imag_abs)


        !
        ! Step 2: handle m=0(radius-local space-time average) using 1D characteristics
        !         based on the perturbation of the local m=0 quantity from the boundary
        !         global space-time average.
        !


        ! Handle radius-local mean variation from global boundary average using 1D characteristics
        call self%compute_boundary_global_average(worker,bc_comm,density_real_m(:,1,1),     &
                                                                 vel1_real_m(:,1,1),        &
                                                                 vel2_real_m(:,1,1),        &
                                                                 vel3_real_m(:,1,1),        &
                                                                 pressure_real_m(:,1,1),    &
                                                                 density_avg,vel1_avg,vel2_avg,vel3_avg,pressure_avg)
        c_avg = sqrt(gam*pressure_avg/density_avg)

        ! Compute perturbation of radius local-avg from boundary-global avg
        ddensity  = density_real_m(:,1,1)  - density_avg 
        dvel1     = vel1_real_m(:,1,1)     - vel1_avg
        dvel2     = vel2_real_m(:,1,1)     - vel2_avg
        dvel3     = vel3_real_m(:,1,1)     - vel3_avg
        dpressure = pressure_real_m(:,1,1) - pressure_avg

        ! Allocate/compute 1d characteristics
        c1_1d = ZERO*ddensity
        c2_1d = ZERO*ddensity
        c3_1d = ZERO*ddensity
        c4_1d = ZERO*ddensity
        c5_1d = ZERO*ddensity
        do iradius = 1,size(ddensity)
            c1_1d(iradius) = -c_avg*c_avg*ddensity(iradius)    +  dpressure(iradius)
            c2_1d(iradius) = density_avg*c_avg*dvel1(iradius)
            c3_1d(iradius) = density_avg*c_avg*dvel2(iradius)
            c4_1d(iradius) = density_avg*c_avg*dvel3(iradius)  +  dpressure(iradius)
            c5_1d(iradius) = -TWO*(pressure_avg - p_user(1))
        end do


        ! Begin composition of radius-local absorbing avg from boundary-global average
        density_real_abs(:,1,1)  = density_avg
        vel1_real_abs(:,1,1)     = vel1_avg
        vel2_real_abs(:,1,1)     = vel2_avg
        vel3_real_abs(:,1,1)     = vel3_avg
        pressure_real_abs(:,1,1) = pressure_avg


        ! Now add local perturbation from the average
        do iradius = 1,size(c1_1d)
            density_real_abs(iradius,1,1)  = density_real_abs(iradius,1,1)  + (-ONE/(c_avg*c_avg))*c1_1d(iradius)  +  (ONE/(TWO*c_avg*c_avg))*c4_1d(iradius)  +  (ONE/(TWO*c_avg*c_avg))*c5_1d(iradius)
            vel1_real_abs(iradius,1,1)     = vel1_real_abs(iradius,1,1)     + (ONE/(density_avg*c_avg))*c2_1d(iradius)
            vel2_real_abs(iradius,1,1)     = vel2_real_abs(iradius,1,1)     + (ONE/(density_avg*c_avg))*c3_1d(iradius)
            vel3_real_abs(iradius,1,1)     = vel3_real_abs(iradius,1,1)     + (ONE/(TWO*density_avg*c_avg))*c4_1d(iradius)  -  (ONE/(TWO*density_avg*c_avg))*c5_1d(iradius)
            pressure_real_abs(iradius,1,1) = pressure_real_abs(iradius,1,1) + HALF*c4_1d(iradius)  +  HALF*c5_1d(iradius)
        end do



    end subroutine compute_absorbing_outlet
    !********************************************************************************






    !>
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine apply_nonreflecting_condition(self,worker,bc_comm, &
                                             density_bar_r, vel1_bar_r, vel2_bar_r, vel3_bar_r, pressure_bar_r, c_bar_r, &
                                             a1_real_m,    a1_imag_m,    &
                                             a2_real_m,    a2_imag_m,    &
                                             a3_real_m,    a3_imag_m,    &
                                             a4_real_m,    a4_imag_m,    &
                                             a5_real_m,    a5_imag_m,    &
                                             a1_real_p,    a1_imag_p,    &
                                             a2_real_p,    a2_imag_p,    &
                                             a3_real_p,    a3_imag_p,    &
                                             a4_real_p,    a4_imag_p,    &
                                             a5_real_p,    a5_imag_p,    &
                                             a1_real,      a1_imag,      &
                                             a2_real,      a2_imag,      &
                                             a3_real,      a3_imag,      &
                                             a4_real,      a4_imag,      &
                                             a5_real,      a5_imag)
        class(outlet_nrbc_lindblad_form3_fix_t),     intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(mpi_comm),             intent(in)      :: bc_comm
        type(AD_D),                 intent(in)      :: density_bar_r(:)
        type(AD_D),                 intent(in)      :: vel1_bar_r(:)
        type(AD_D),                 intent(in)      :: vel2_bar_r(:)
        type(AD_D),                 intent(in)      :: vel3_bar_r(:)
        type(AD_D),                 intent(in)      :: pressure_bar_r(:)
        type(AD_D),                 intent(in)      :: c_bar_r(:)
        type(AD_D),                 intent(inout)   :: a1_real_m(:,:,:)
        type(AD_D),                 intent(inout)   :: a1_imag_m(:,:,:)
        type(AD_D),                 intent(inout)   :: a2_real_m(:,:,:)
        type(AD_D),                 intent(inout)   :: a2_imag_m(:,:,:)
        type(AD_D),                 intent(inout)   :: a3_real_m(:,:,:)
        type(AD_D),                 intent(inout)   :: a3_imag_m(:,:,:)
        type(AD_D),                 intent(inout)   :: a4_real_m(:,:,:)
        type(AD_D),                 intent(inout)   :: a4_imag_m(:,:,:)
        type(AD_D),                 intent(inout)   :: a5_real_m(:,:,:)
        type(AD_D),                 intent(inout)   :: a5_imag_m(:,:,:)
        type(AD_D),                 intent(inout)   :: a1_real_p(:,:,:)
        type(AD_D),                 intent(inout)   :: a1_imag_p(:,:,:)
        type(AD_D),                 intent(inout)   :: a2_real_p(:,:,:)
        type(AD_D),                 intent(inout)   :: a2_imag_p(:,:,:)
        type(AD_D),                 intent(inout)   :: a3_real_p(:,:,:)
        type(AD_D),                 intent(inout)   :: a3_imag_p(:,:,:)
        type(AD_D),                 intent(inout)   :: a4_real_p(:,:,:)
        type(AD_D),                 intent(inout)   :: a4_imag_p(:,:,:)
        type(AD_D),                 intent(inout)   :: a5_real_p(:,:,:)
        type(AD_D),                 intent(inout)   :: a5_imag_p(:,:,:)
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
        
        integer(ik) :: iradius, itheta, ntheta, itime

        type(AD_D)  :: beta, B3_real, B3_imag, B4_real, B4_imag, &
                       density_bar, vel1_bar, vel2_bar, vel3_bar, pressure_bar, c_bar


        ! Outoing amplitudes from interior
        a1_real = a1_real_m
        a1_imag = a1_imag_m
        a2_real = a2_real_m
        a2_imag = a2_imag_m
        a3_real = a3_real_m
        a3_imag = a3_imag_m
        a4_real = a4_real_m
        a4_imag = a4_imag_m

        ! Incoming amplitudes from exterior
        a5_real = ZERO*a5_real_p
        a5_imag = ZERO*a5_imag_p


!        itime = 1   ! Time-constant
!        do iradius = 1,size(a1_real_m,1)
!            ! Get average parts
!            density_bar  = density_bar_r(iradius)
!            vel1_bar     = vel1_bar_r(iradius)
!            vel2_bar     = vel2_bar_r(iradius)
!            vel3_bar     = vel3_bar_r(iradius)
!            pressure_bar = pressure_bar_r(iradius)
!            c_bar        = sqrt(gam*pressure_bar/density_bar)
!
!            ! starting with 2 here because the first mode is treated with 1D characteristics
!            ntheta = size(a1_real_m,2)
!            do itheta = 2,ntheta
!                ! Account for sign(mode) in the calculation of beta. The second half of the
!                ! modes are negative frequencies.
!                if (itheta <= (ntheta-1)/2 + 1) then
!                    beta = sqrt(c_bar*c_bar  -  (vel3_bar*vel3_bar + vel2_bar*vel2_bar))
!                else if (itheta > (ntheta-1)/2 + 1) then
!                    beta = -sqrt(c_bar*c_bar  -  (vel3_bar*vel3_bar + vel2_bar*vel2_bar))
!                end if
!
!                ! The imaginary part of beta has already been accounted for in
!                ! the expressions for B2 and B3
!                B3_real = -TWO*vel3_bar*vel2_bar/(vel2_bar*vel2_bar + beta*beta)
!                B3_imag = -TWO*beta*vel3_bar/(vel2_bar*vel2_bar + beta*beta)
!
!                B4_real = (beta*beta - vel2_bar*vel2_bar)/(beta*beta + vel2_bar*vel2_bar)
!                B4_imag = -TWO*beta*vel2_bar/(beta*beta + vel2_bar*vel2_bar)
!
!                a5_real(iradius,itheta,itime) = (B3_real*a3_real_m(iradius,itheta,itime) - B3_imag*a3_imag_m(iradius,itheta,itime))  &   ! A3*c3 (real)
!                                              - (B4_real*a4_real_m(iradius,itheta,itime) - B4_imag*a4_imag_m(iradius,itheta,itime))      ! A4*c4 (real)
!                a5_imag(iradius,itheta,itime) = (B3_imag*a3_real_m(iradius,itheta,itime) + B3_real*a3_imag_m(iradius,itheta,itime))  &   ! A3*c3 (imag)
!                                              - (B4_imag*a4_real_m(iradius,itheta,itime) + B4_real*a4_imag_m(iradius,itheta,itime))      ! A4*c4 (imag)
!            end do !itheta
!        end do !iradius




    end subroutine apply_nonreflecting_condition
    !********************************************************************************



    !>  Compute boundary global average by averaging spatio-temporal average over
    !!  radial stations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/16/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_boundary_global_average(self,worker,bc_comm,density_bar,vel1_bar,vel2_bar,vel3_bar,pressure_bar, &
                                                                   density_avg,vel1_avg,vel2_avg,vel3_avg,pressure_avg)
        class(outlet_nrbc_lindblad_form3_fix_t),  intent(inout)   :: self
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
                density_avg  = density_avg  + dr*(density_bar(irad+1)  + density_bar(irad))/TWO
                vel1_avg     = vel1_avg     + dr*(vel1_bar(irad+1)     + vel1_bar(irad))/TWO
                vel2_avg     = vel2_avg     + dr*(vel2_bar(irad+1)     + vel2_bar(irad))/TWO
                vel3_avg     = vel3_avg     + dr*(vel3_bar(irad+1)     + vel3_bar(irad))/TWO
                pressure_avg = pressure_avg + dr*(pressure_bar(irad+1) + pressure_bar(irad))/TWO
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

    end subroutine compute_boundary_global_average
    !********************************************************************************







end module bc_state_outlet_nrbc_lindblad_form3_fix
