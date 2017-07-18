module bc_state_wall
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO, ONE
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use mpi_f08,                only: mpi_comm
    use DNAD_D
    implicit none
    


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: wall_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type wall_t
    !*******************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(wall_t),   intent(inout) :: self
        

        !
        ! Set operator name
        !
        call self%set_name('Wall')
        call self%set_family('Wall')


!        !
!        ! Set operator equations
!        !
!        call self%set_equation('Density'   )
!        call self%set_equation('Momentum-1')
!        call self%set_equation('Momentum-2')
!        call self%set_equation('Momentum-3')
!        call self%set_equation('Energy'    )


    end subroutine init
    !********************************************************************************














    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(wall_t),          intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        class(properties_t),    intent(inout)   :: prop
        type(mpi_comm),         intent(in)      :: bc_COMM


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                      &
            density_m,  mom1_m,  mom2_m,  mom3_m,  energy_m,            &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc,           &
            drho_dx_m, drhou_dx_m, drhov_dx_m, drhow_dx_m, drhoE_dx_m,  &
            drho_dy_m, drhou_dy_m, drhov_dy_m, drhow_dy_m, drhoE_dy_m,  &
            drho_dz_m, drhou_dz_m, drhov_dz_m, drhow_dz_m, drhoE_dz_m,  &
            u_m, v_m, w_m, p_m, normal_velocity, normal_velocity_1, normal_velocity_2, normal_velocity_3

        real(rk),   allocatable, dimension(:)   :: unorm_1, unorm_2, unorm_3, r
    

        real(rk), allocatable, dimension(:) ::      &
           u_grid, v_grid, w_grid, det_jacobian_grid

        type(AD_D), allocatable, dimension(:,:) :: grad_density, grad_mom1, grad_mom2, grad_mom3, grad_energy

        det_jacobian_grid = worker%get_det_jacobian_grid_face('value')


        u_grid = worker%get_grid_velocity_face('u_grid')
        v_grid = worker%get_grid_velocity_face('v_grid')
        w_grid = worker%get_grid_velocity_face('w_grid')

        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_m = worker%get_primary_field_face('Density'   , 'value', 'face interior')
        mom1_m    = worker%get_primary_field_face('Momentum-1', 'value', 'face interior')
        mom2_m    = worker%get_primary_field_face('Momentum-2', 'value', 'face interior')
        mom3_m    = worker%get_primary_field_face('Momentum-3', 'value', 'face interior')
        energy_m  = worker%get_primary_field_face('Energy'    , 'value', 'face interior')

        density_m = density_m/det_jacobian_grid
        mom1_m = mom1_m/det_jacobian_grid
        mom2_m = mom2_m/det_jacobian_grid
        mom3_m = mom3_m/det_jacobian_grid
        energy_m = energy_m/det_jacobian_grid



        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        r = worker%coordinate('1','boundary')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / r
        else if (worker%coordinate_system() == 'Cartesian') then

        else
            call chidg_signal(FATAL,"inlet, bad coordinate system")
        end if
    
        !
        ! BUG!!! need to be getting gradient without lift b/c it is not computed at boundaries!!!!!
        !
        grad_density    = worker%get_primary_field_grad_ale_face('Density'   , 'gradient', 'face interior')
        grad_mom1       = worker%get_primary_field_grad_ale_face('Momentum-1', 'gradient', 'face interior')
        grad_mom2       = worker%get_primary_field_grad_ale_face('Momentum-2', 'gradient', 'face interior')
        grad_mom3       = worker%get_primary_field_grad_ale_face('Momentum-3', 'gradient', 'face interior')
        grad_energy     = worker%get_primary_field_grad_ale_face('Energy'    , 'gradient', 'face interior')



        drho_dx_m  = grad_density(:,1) !worker%get_primary_field_face('Density'   , 'grad1', 'face interior')
        drho_dy_m  = grad_density(:,2) !worker%get_primary_field_face('Density'   , 'grad2', 'face interior')
        drho_dz_m  = grad_density(:,3) !worker%get_primary_field_face('Density'   , 'grad3', 'face interior')

        drhou_dx_m = grad_mom1(:,1)!worker%get_primary_field_face('Momentum-1', 'grad1', 'face interior')
        drhou_dy_m = grad_mom1(:,2)!worker%get_primary_field_face('Momentum-1', 'grad2', 'face interior')
        drhou_dz_m = grad_mom1(:,3)!worker%get_primary_field_face('Momentum-1', 'grad3', 'face interior')

        drhov_dx_m = grad_mom2(:,1)!worker%get_primary_field_face('Momentum-2', 'grad1', 'face interior')
        drhov_dy_m = grad_mom2(:,2)!worker%get_primary_field_face('Momentum-2', 'grad2', 'face interior')
        drhov_dz_m = grad_mom2(:,3)!worker%get_primary_field_face('Momentum-2', 'grad3', 'face interior')

        drhow_dx_m = grad_mom3(:,1)!worker%get_primary_field_face('Momentum-3', 'grad1', 'face interior')
        drhow_dy_m = grad_mom3(:,2)!worker%get_primary_field_face('Momentum-3', 'grad2', 'face interior')
        drhow_dz_m = grad_mom3(:,3)!worker%get_primary_field_face('Momentum-3', 'grad3', 'face interior')
        
        drhoE_dx_m = grad_energy(:,1)!worker%get_primary_field_face('Energy'    , 'grad1', 'face interior')
        drhoE_dy_m = grad_energy(:,2)!worker%get_primary_field_face('Energy'    , 'grad2', 'face interior')
        drhoE_dz_m = grad_energy(:,3)!worker%get_primary_field_face('Energy'    , 'grad3', 'face interior')



!        drho_dx_m  = worker%get_primary_field_face('Density'   , 'grad1', 'face interior')
!        drho_dy_m  = worker%get_primary_field_face('Density'   , 'grad2', 'face interior')
!        drho_dz_m  = worker%get_primary_field_face('Density'   , 'grad3', 'face interior')
!
!        drhou_dx_m = worker%get_primary_field_face('Momentum-1', 'grad1', 'face interior')
!        drhou_dy_m = worker%get_primary_field_face('Momentum-1', 'grad2', 'face interior')
!        drhou_dz_m = worker%get_primary_field_face('Momentum-1', 'grad3', 'face interior')
!
!        drhov_dx_m = worker%get_primary_field_face('Momentum-2', 'grad1', 'face interior')
!        drhov_dy_m = worker%get_primary_field_face('Momentum-2', 'grad2', 'face interior')
!        drhov_dz_m = worker%get_primary_field_face('Momentum-2', 'grad3', 'face interior')
!
!        drhow_dx_m = worker%get_primary_field_face('Momentum-3', 'grad1', 'face interior')
!        drhow_dy_m = worker%get_primary_field_face('Momentum-3', 'grad2', 'face interior')
!        drhow_dz_m = worker%get_primary_field_face('Momentum-3', 'grad3', 'face interior')
!        
!        drhoE_dx_m = worker%get_primary_field_face('Energy'    , 'grad1', 'face interior')
!        drhoE_dy_m = worker%get_primary_field_face('Energy'    , 'grad2', 'face interior')
!        drhoE_dz_m = worker%get_primary_field_face('Energy'    , 'grad3', 'face interior')


        !
        ! Get normal vectors
        !
        unorm_1 = worker%unit_normal(1)
        unorm_2 = worker%unit_normal(2)
        unorm_3 = worker%unit_normal(3)



        ! Initialize arrays
        density_bc = density_m
        mom1_bc    = density_m
        mom2_bc    = density_m
        mom3_bc    = density_m
        energy_bc  = density_m


        ! Zero momentum
        !mom1_bc = ZERO
        !mom2_bc = ZERO
        !mom3_bc = ZERO


        ! Set relative velocity to zero, ie
        ! set fluid velocity equal to grid/wall velocity.
        mom1_bc = density_m*u_grid
        mom2_bc = density_m*v_grid
        mom3_bc = density_m*w_grid
        !
        ! We want:  W dot n = 0
        !
        u_m = mom1_m/density_m-u_grid
        v_m = mom2_m/density_m-v_grid
        w_m = mom3_m/density_m-w_grid


        !
        ! Energy subtract change in kinetic energy
        !
        energy_bc = energy_m - (density_m*HALF)*(u_m*u_m  +  v_m*v_m  +  w_m*w_m)



        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_bc = mom2_bc * r
        else if (worker%coordinate_system() == 'Cartesian') then

        else
            call chidg_signal(FATAL,"inlet, bad coordinate system")
        end if



        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density'   , density_bc, 'value')
        call worker%store_bc_state('Momentum-1', mom1_bc,    'value')
        call worker%store_bc_state('Momentum-2', mom2_bc,    'value')
        call worker%store_bc_state('Momentum-3', mom3_bc,    'value')
        call worker%store_bc_state('Energy'    , energy_bc,  'value')


        drho_dx_m = ZERO
        drho_dy_m = ZERO
        drho_dz_m = ZERO
        call worker%store_bc_state('Density'   , drho_dx_m, 'grad1')
        call worker%store_bc_state('Density'   , drho_dy_m, 'grad2')
        call worker%store_bc_state('Density'   , drho_dz_m, 'grad3')
                                                
        call worker%store_bc_state('Momentum-1', drhou_dx_m, 'grad1')
        call worker%store_bc_state('Momentum-1', drhou_dy_m, 'grad2')
        call worker%store_bc_state('Momentum-1', drhou_dz_m, 'grad3')
                                                
        call worker%store_bc_state('Momentum-2', drhov_dx_m, 'grad1')
        call worker%store_bc_state('Momentum-2', drhov_dy_m, 'grad2')
        call worker%store_bc_state('Momentum-2', drhov_dz_m, 'grad3')
                                                
        call worker%store_bc_state('Momentum-3', drhow_dx_m, 'grad1')
        call worker%store_bc_state('Momentum-3', drhow_dy_m, 'grad2')
        call worker%store_bc_state('Momentum-3', drhow_dz_m, 'grad3')

        drhoE_dx_m = ZERO
        drhoE_dy_m = ZERO
        drhoE_dz_m = ZERO
        call worker%store_bc_state('Energy'    , drhoE_dx_m, 'grad1')
        call worker%store_bc_state('Energy'    , drhoE_dy_m, 'grad2')
        call worker%store_bc_state('Energy'    , drhoE_dz_m, 'grad3')



    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_wall
