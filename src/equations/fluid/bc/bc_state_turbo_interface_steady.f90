module bc_state_turbo_interface_steady
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF, ME, CYLINDRICAL,    &
                                      XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, PI
    use mod_fluid,              only: gam, Rgas, cp

    use type_point,             only: point_t
    use type_mesh,              only: mesh_t
    use type_bc_state,          only: bc_state_t
    use bc_giles_HB_base,       only: giles_HB_base_t
    use type_bc_patch,          only: bc_patch_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_face_info,         only: face_info_t, face_info_constructor
    use type_element_info,      only: element_info_t
    use ieee_arithmetic,        only: ieee_is_nan
    use mpi_f08,                only: mpi_comm
    use DNAD_D
    implicit none



    !>  Name: inlet - 3D Giles
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
    type, public, extends(giles_HB_base_t) :: turbo_interface_steady_t

    contains

        procedure   :: init                         ! Set-up bc state with options/name etc.
        procedure   :: compute_bc_state             ! boundary condition function implementation
        procedure   :: compute_absorbing_interface
        procedure   :: apply_nonreflecting_condition

    end type turbo_interface_steady_t
    !*********************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. average_pressure 
    !!  @date   2/8/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(turbo_interface_steady_t),   intent(inout) :: self

        ! Set name, family
        call self%set_name('Steady Turbo Interface')
        call self%set_family('Inlet')

        ! Add functions
        call self%bcproperties%add('Pitch A', 'Required')
        call self%bcproperties%add('Pitch B', 'Required')

!        call self%bcproperties%add('Spatial Periodicity A', 'Required')
!        call self%bcproperties%add('Spatial Periodicity A', 'Required')

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
        class(turbo_interface_steady_t),   intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop
        type(mpi_comm),                             intent(in)      :: bc_comm


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                                                      &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc, pressure_bc, vel1_bc, vel2_bc, vel3_bc,   &
            grad1_density_m, grad1_mom1_m, grad1_mom2_m, grad1_mom3_m, grad1_energy_m,                  &
            grad2_density_m, grad2_mom1_m, grad2_mom2_m, grad2_mom3_m, grad2_energy_m,                  &
            grad3_density_m, grad3_mom1_m, grad3_mom2_m, grad3_mom3_m, grad3_energy_m

        type(AD_D), allocatable, dimension(:,:) ::                                                                      &
            density_check_real_gq, vel1_check_real_gq, vel2_check_real_gq, vel3_check_real_gq, pressure_check_real_gq,  &
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
            density_hat_real_gq, vel1_hat_real_gq, vel2_hat_real_gq, vel3_hat_real_gq, pressure_hat_real_gq,                        &
            density_hat_imag_gq, vel1_hat_imag_gq, vel2_hat_imag_gq, vel3_hat_imag_gq, pressure_hat_imag_gq,                        &
            density_grid_m, vel1_grid_m, vel2_grid_m, vel3_grid_m, pressure_grid_m, c_grid_m,                                       &
            density_grid_p, vel1_grid_p, vel2_grid_p, vel3_grid_p, pressure_grid_p, c_grid_p

        real(rk), allocatable, dimension(:)   :: r
        character(1)    :: side

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
        call self%compute_spatial_dft(worker,bc_comm,'A',                            &
                                      density_check_real_m,  density_check_imag_m,   &
                                      vel1_check_real_m,     vel1_check_imag_m,      &
                                      vel2_check_real_m,     vel2_check_imag_m,      &
                                      vel3_check_real_m,     vel3_check_imag_m,      &
                                      pressure_check_real_m, pressure_check_imag_m,  &
                                      c_check_real_m,        c_check_imag_m,         &
                                      density_hat_real_m,    density_hat_imag_m,     &
                                      vel1_hat_real_m,       vel1_hat_imag_m,        &
                                      vel2_hat_real_m,       vel2_hat_imag_m,        &
                                      vel3_hat_real_m,       vel3_hat_imag_m,        &
                                      pressure_hat_real_m,   pressure_hat_imag_m,    &
                                      c_hat_real_m,          c_hat_imag_m)





        ! Get exterior perturbation
        call self%get_q_side(worker,bc_comm,'B', &
                             density_grid_p,     &
                             vel1_grid_p,        &
                             vel2_grid_p,        &
                             vel3_grid_p,        &
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
        call self%compute_spatial_dft(worker,bc_comm,'B',                               &
                                      density_check_real_p,     density_check_imag_p,   &
                                      vel1_check_real_p,        vel1_check_imag_p,      &
                                      vel2_check_real_p,        vel2_check_imag_p,      &
                                      vel3_check_real_p,        vel3_check_imag_p,      &
                                      pressure_check_real_p,    pressure_check_imag_p,  &
                                      c_check_real_p,           c_check_imag_p,         &
                                      density_hat_real_p,       density_hat_imag_p,     &
                                      vel1_hat_real_p,          vel1_hat_imag_p,        &
                                      vel2_hat_real_p,          vel2_hat_imag_p,        &
                                      vel3_hat_real_p,          vel3_hat_imag_p,        &
                                      pressure_hat_real_p,      pressure_hat_imag_p,    &
                                      c_hat_real_p,             c_hat_imag_p)



        ! Compute q_abs = f(q_p,q_m)
        side = self%get_face_side(worker)
        call self%compute_absorbing_interface(worker,bc_comm,side,                          &
                                              density_hat_real_m,    density_hat_imag_m,    &
                                              vel1_hat_real_m,       vel1_hat_imag_m,       &
                                              vel2_hat_real_m,       vel2_hat_imag_m,       &
                                              vel3_hat_real_m,       vel3_hat_imag_m,       &
                                              pressure_hat_real_m,   pressure_hat_imag_m,   &
                                              c_hat_real_m,          c_hat_imag_m,          &
                                              density_hat_real_p,    density_hat_imag_p,    &
                                              vel1_hat_real_p,       vel1_hat_imag_p,       &
                                              vel2_hat_real_p,       vel2_hat_imag_p,       &
                                              vel3_hat_real_p,       vel3_hat_imag_p,       &
                                              pressure_hat_real_p,   pressure_hat_imag_p,   &
                                              c_hat_real_p,          c_hat_imag_p,          &
                                              density_hat_real_abs,  density_hat_imag_abs,  &
                                              vel1_hat_real_abs,     vel1_hat_imag_abs,     &
                                              vel2_hat_real_abs,     vel2_hat_imag_abs,     &
                                              vel3_hat_real_abs,     vel3_hat_imag_abs,     &
                                              pressure_hat_real_abs, pressure_hat_imag_abs)


        ! q_abs(r_gq) = I(q_abs(r_aux))
        call self%interpolate_raux_to_rgq(worker,bc_comm,   &
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
        call self%compute_spatial_idft_gq(worker,bc_comm,side,                              &
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
        ! Store boundary condition state. q_bc
        !
        call worker%store_bc_state('Density'   , density_bc, 'value')
        call worker%store_bc_state('Momentum-1', mom1_bc,    'value')
        call worker%store_bc_state('Momentum-2', mom2_bc,    'value')
        call worker%store_bc_state('Momentum-3', mom3_bc,    'value')
        call worker%store_bc_state('Energy'    , energy_bc,  'value')



    end subroutine compute_bc_state
    !*********************************************************************************





    !> DANIEL'S FORMULATION
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    subroutine compute_absorbing_interface(self,worker,bc_comm,side,            &
                                           density_real_m,    density_imag_m,   &
                                           vel1_real_m,       vel1_imag_m,      &
                                           vel2_real_m,       vel2_imag_m,      &
                                           vel3_real_m,       vel3_imag_m,      &
                                           pressure_real_m,   pressure_imag_m,  &
                                           c_real_m,          c_imag_m,         &
                                           density_real_p,    density_imag_p,   &
                                           vel1_real_p,       vel1_imag_p,      &
                                           vel2_real_p,       vel2_imag_p,      &
                                           vel3_real_p,       vel3_imag_p,      &
                                           pressure_real_p,   pressure_imag_p,  &
                                           c_real_p,          c_imag_p,         &
                                           density_real_abs,  density_imag_abs, &
                                           vel1_real_abs,     vel1_imag_abs,    &
                                           vel2_real_abs,     vel2_imag_abs,    &
                                           vel3_real_abs,     vel3_imag_abs,    &
                                           pressure_real_abs, pressure_imag_abs)
        class(turbo_interface_steady_t),    intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        type(mpi_comm),                     intent(in)      :: bc_comm
        character(1),                       intent(in)      :: side
        type(AD_D),                         intent(inout)   :: density_real_m(:,:,:)
        type(AD_D),                         intent(inout)   :: density_imag_m(:,:,:)
        type(AD_D),                         intent(inout)   :: vel1_real_m(:,:,:)
        type(AD_D),                         intent(inout)   :: vel1_imag_m(:,:,:)
        type(AD_D),                         intent(inout)   :: vel2_real_m(:,:,:)
        type(AD_D),                         intent(inout)   :: vel2_imag_m(:,:,:)
        type(AD_D),                         intent(inout)   :: vel3_real_m(:,:,:)
        type(AD_D),                         intent(inout)   :: vel3_imag_m(:,:,:)
        type(AD_D),                         intent(inout)   :: pressure_real_m(:,:,:)
        type(AD_D),                         intent(inout)   :: pressure_imag_m(:,:,:)
        type(AD_D),                         intent(inout)   :: c_real_m(:,:,:)
        type(AD_D),                         intent(inout)   :: c_imag_m(:,:,:)
        type(AD_D),                         intent(inout)   :: density_real_p(:,:,:)
        type(AD_D),                         intent(inout)   :: density_imag_p(:,:,:)
        type(AD_D),                         intent(inout)   :: vel1_real_p(:,:,:)
        type(AD_D),                         intent(inout)   :: vel1_imag_p(:,:,:)
        type(AD_D),                         intent(inout)   :: vel2_real_p(:,:,:)
        type(AD_D),                         intent(inout)   :: vel2_imag_p(:,:,:)
        type(AD_D),                         intent(inout)   :: vel3_real_p(:,:,:)
        type(AD_D),                         intent(inout)   :: vel3_imag_p(:,:,:)
        type(AD_D),                         intent(inout)   :: pressure_real_p(:,:,:)
        type(AD_D),                         intent(inout)   :: pressure_imag_p(:,:,:)
        type(AD_D),                         intent(inout)   :: c_real_p(:,:,:)
        type(AD_D),                         intent(inout)   :: c_imag_p(:,:,:)
        type(AD_D), allocatable,            intent(inout)   :: density_real_abs(:,:,:)
        type(AD_D), allocatable,            intent(inout)   :: density_imag_abs(:,:,:)
        type(AD_D), allocatable,            intent(inout)   :: vel1_real_abs(:,:,:)
        type(AD_D), allocatable,            intent(inout)   :: vel1_imag_abs(:,:,:)
        type(AD_D), allocatable,            intent(inout)   :: vel2_real_abs(:,:,:)
        type(AD_D), allocatable,            intent(inout)   :: vel2_imag_abs(:,:,:)
        type(AD_D), allocatable,            intent(inout)   :: vel3_real_abs(:,:,:)
        type(AD_D), allocatable,            intent(inout)   :: vel3_imag_abs(:,:,:)
        type(AD_D), allocatable,            intent(inout)   :: pressure_real_abs(:,:,:)
        type(AD_D), allocatable,            intent(inout)   :: pressure_imag_abs(:,:,:)

        type(AD_D), allocatable, dimension(:,:,:) ::                &
            a1_real, a2_real, a3_real, a4_real, a5_real,            &
            a1_imag, a2_imag, a3_imag, a4_imag, a5_imag,            &
            a1_real_m, a2_real_m, a3_real_m, a4_real_m, a5_real_m,  &
            a1_imag_m, a2_imag_m, a3_imag_m, a4_imag_m, a5_imag_m,  &
            a1_real_p, a2_real_p, a3_real_p, a4_real_p, a5_real_p,  &
            a1_imag_p, a2_imag_p, a3_imag_p, a4_imag_p, a5_imag_p

        real(rk),   allocatable, dimension(:,:) :: unorm

        print*, 'WARNING: Inconsistent use of Pitch A in eigenvalue calc'

        ! Project interior to eigenmodes
        call self%primitive_to_eigenmodes(worker,bc_comm,                   &
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
        call self%primitive_to_eigenmodes(worker,bc_comm,                   &
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


        ! Allocate and zero absorbing modes
        a1_real = ZERO*a1_real_m
        a1_imag = ZERO*a1_real_m
        a2_real = ZERO*a1_real_m
        a2_imag = ZERO*a1_real_m
        a3_real = ZERO*a1_real_m
        a3_imag = ZERO*a1_real_m
        a4_real = ZERO*a1_real_m
        a4_imag = ZERO*a1_real_m
        a5_real = ZERO*a1_real_m
        a5_imag = ZERO*a1_real_m

        call self%apply_nonreflecting_condition(worker,bc_comm,side,        &
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



        ! To initialize average and storage
        density_real_abs  = density_real_m
        density_imag_abs  = density_imag_m
        vel1_real_abs     = vel1_real_m
        vel1_imag_abs     = vel1_imag_m
        vel2_real_abs     = vel2_real_m
        vel2_imag_abs     = vel2_imag_m
        vel3_real_abs     = vel3_real_m
        vel3_imag_abs     = vel3_imag_m
        pressure_real_abs = pressure_real_m
        pressure_imag_abs = pressure_imag_m


        ! Convert back to primitive variables
        call self%eigenmodes_to_primitive(worker,bc_comm,                       &
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


        ! Update space-time average as average of 'A' and 'B'
        density_real_abs(:,1,1)  = (density_real_m(:,1,1)  + density_real_p(:,1,1))/TWO
        vel1_real_abs(:,1,1)     = (vel1_real_m(:,1,1)     + vel1_real_p(:,1,1))/TWO
        vel2_real_abs(:,1,1)     = (vel2_real_m(:,1,1)     + vel2_real_p(:,1,1))/TWO
        vel3_real_abs(:,1,1)     = (vel3_real_m(:,1,1)     + vel3_real_p(:,1,1))/TWO
        pressure_real_abs(:,1,1) = (pressure_real_m(:,1,1) + pressure_real_p(:,1,1))/TWO

        ! Zero imaginary part
        density_imag_abs(:,1,1)  = ZERO
        vel1_imag_abs(:,1,1)     = ZERO
        vel2_imag_abs(:,1,1)     = ZERO
        vel3_imag_abs(:,1,1)     = ZERO
        pressure_imag_abs(:,1,1) = ZERO

                                       
    end subroutine compute_absorbing_interface
    !********************************************************************************





    !>
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine apply_nonreflecting_condition(self,worker,bc_comm,side, &
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
        class(turbo_interface_steady_t),    intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        type(mpi_comm),             intent(in)      :: bc_comm
        character(1),               intent(in)      :: side
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


        ! NOTE the amplitude strategy here is for 1d characteristics. The unsteady amplitudes
        ! have diffferent propagation directions so the logic will be different. 


        if (side=='A') then

            ! Extrapolate amplitudes from upstream
            a1_real = a1_real_m
            a1_imag = a1_imag_m
            a2_real = a2_real_m
            a2_imag = a2_imag_m
            a3_real = a3_real_m
            a3_imag = a3_imag_m
            a4_real = a4_real_m
            a4_imag = a4_imag_m

            ! Zero amplitudes from downstream
            a5_real = ZERO*a5_real_p
            a5_imag = ZERO*a5_imag_p

            ! Adjust steady modes based on Giles' original steady formulation.
            itime = 1   ! Time-constant
            do iradius = 1,size(a1_real_m,1)
                ! Get average parts
                density_bar  = density_bar_r(iradius)
                vel1_bar     = vel1_bar_r(iradius)
                vel2_bar     = vel2_bar_r(iradius)
                vel3_bar     = vel3_bar_r(iradius)
                pressure_bar = pressure_bar_r(iradius)
                c_bar        = sqrt(gam*pressure_bar/density_bar)

                ! starting with 2 here because the first mode is treated with 1D characteristics
                ntheta = size(a1_real_m,2)
                do itheta = 2,ntheta
                    ! Account for sign(mode) in the calculation of beta. The second half of the
                    ! modes are negative frequencies.
                    if (itheta <= (ntheta-1)/2 + 1) then
                        beta = sqrt(c_bar*c_bar  -  (vel3_bar*vel3_bar + vel2_bar*vel2_bar))
                    else if (itheta > (ntheta-1)/2 + 1) then
                        beta = -sqrt(c_bar*c_bar  -  (vel3_bar*vel3_bar + vel2_bar*vel2_bar))
                    end if

                    ! The imaginary part of beta has already been accounted for in
                    ! the expressions for B2 and B3
                    B3_real = -TWO*vel3_bar*vel2_bar/(vel2_bar*vel2_bar + beta*beta)
                    B3_imag = -TWO*beta*vel3_bar/(vel2_bar*vel2_bar + beta*beta)

                    B4_real = (beta*beta - vel2_bar*vel2_bar)/(beta*beta + vel2_bar*vel2_bar)
                    B4_imag = -TWO*beta*vel2_bar/(beta*beta + vel2_bar*vel2_bar)

                    a5_real(iradius,itheta,itime) = (B3_real*a3_real_m(iradius,itheta,itime) - B3_imag*a3_imag_m(iradius,itheta,itime))  &   ! A3*c3 (real)
                                                  - (B4_real*a4_real_m(iradius,itheta,itime) - B4_imag*a4_imag_m(iradius,itheta,itime))      ! A4*c4 (real)
                    a5_imag(iradius,itheta,itime) = (B3_imag*a3_real_m(iradius,itheta,itime) + B3_real*a3_imag_m(iradius,itheta,itime))  &   ! A3*c3 (imag)
                                                  - (B4_imag*a4_real_m(iradius,itheta,itime) + B4_real*a4_imag_m(iradius,itheta,itime))      ! A4*c4 (imag)
                end do !itheta
            end do !iradius


        else if (side=='B') then


            ! Zero amplitudes from upsteady
!            a1_real(:,:,1) = ZERO
!            a1_imag(:,:,1) = ZERO
!            a2_real(:,:,1) = ZERO
!            a2_imag(:,:,1) = ZERO
!            a3_real(:,:,1) = ZERO
!            a3_imag(:,:,1) = ZERO
!            a4_real(:,:,1) = ZERO
!            a4_imag(:,:,1) = ZERO
!
!            a5_real(:,:,1) = a5_real_p(:,:,1)
!            a5_imag(:,:,1) = a5_imag_p(:,:,1)

            ! Zero incoming amplitudes
            a1_real = ZERO
            a1_imag = ZERO
            a2_real = ZERO
            a2_imag = ZERO
            a3_real = ZERO
            a3_imag = ZERO
            a4_real = ZERO
            a4_imag = ZERO

            ! Extrapolate outgoing amplitudes
            a5_real = a5_real_p
            a5_imag = a5_imag_p


            itime = 1   ! Time-constant
            do iradius = 1,size(a1_real_m,1)
                ! Get average parts
                density_bar  = density_bar_r(iradius)
                vel1_bar     = vel1_bar_r(iradius)
                vel2_bar     = vel2_bar_r(iradius)
                vel3_bar     = vel3_bar_r(iradius)
                pressure_bar = pressure_bar_r(iradius)
                c_bar        = sqrt(gam*pressure_bar/density_bar)

                ! starting with 2 here because the first mode is treated with 1D characteristics
                ntheta = size(a1_real_m,2)
                do itheta = 2,ntheta
                    ! Account for sign(mode) in the calculation of beta. The second half of the
                    ! modes are negative frequencies.
                    if (itheta <= (ntheta-1)/2 + 1) then
                        beta = sqrt(c_bar*c_bar  -  (vel3_bar*vel3_bar + vel2_bar*vel2_bar))
                    else if (itheta > (ntheta-1)/2 + 1) then
                        beta = -sqrt(c_bar*c_bar  -  (vel3_bar*vel3_bar + vel2_bar*vel2_bar))
                    end if


                    B3_real = -vel2_bar/(c_bar + vel3_bar)
                    B3_imag = -beta/(c_bar + vel3_bar)

                    B4_real = (vel2_bar*vel2_bar - beta*beta)/((c_bar + vel3_bar)**TWO)
                    B4_imag = TWO*vel2_bar*beta/((c_bar + vel3_bar)**TWO)


                    a1_real(iradius,itheta,itime) = ZERO
                    a1_imag(iradius,itheta,itime) = ZERO
                    a2_real(iradius,itheta,itime) = ZERO
                    a2_imag(iradius,itheta,itime) = ZERO
                                          
                    a3_real(iradius,itheta,itime) = (B3_real*a5_real_p(iradius,itheta,itime) - B3_imag*a5_imag_p(iradius,itheta,itime))
                    a3_imag(iradius,itheta,itime) = (B3_imag*a5_real_p(iradius,itheta,itime) + B3_real*a5_imag_p(iradius,itheta,itime))
                                                                                                                 
                    a4_real(iradius,itheta,itime) = (B4_real*a5_real_p(iradius,itheta,itime) - B4_imag*a5_imag_p(iradius,itheta,itime))
                    a4_imag(iradius,itheta,itime) = (B4_imag*a5_real_p(iradius,itheta,itime) + B4_real*a5_imag_p(iradius,itheta,itime))

                end do !itheta
            end do !iradius

        else
            call chidg_signal(FATAL,"turbo_interface_steady%apply_nonreflecting_condition: invalid input for argument 'side': 'A' or 'B'")
        end if



    end subroutine apply_nonreflecting_condition
    !********************************************************************************






















end module bc_state_turbo_interface_steady
