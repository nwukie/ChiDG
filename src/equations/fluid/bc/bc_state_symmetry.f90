module bc_state_symmetry
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO
    use mod_fluid,              only: omega
    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use mpi_f08,                only: mpi_comm
    use ieee_arithmetic
    use DNAD_D
    implicit none
    


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: symmetry_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type symmetry_t
    !*******************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(symmetry_t),   intent(inout) :: self
        
        ! Set operator name
        call self%set_name("Symmetry")
        call self%set_family("Symmetry")


    end subroutine init
    !********************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(symmetry_t),      intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        class(properties_t),    intent(inout)   :: prop
        type(mpi_comm),         intent(in)      :: bc_COMM

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                      &
            density_m,  mom1_m,  mom2_m,  mom3_m,  energy_m,            &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc,           &
            grad1_density_m, grad1_mom1_m, grad1_mom2_m, grad1_mom3_m, grad1_energy_m,  &
            grad2_density_m, grad2_mom1_m, grad2_mom2_m, grad2_mom3_m, grad2_energy_m,  &
            grad3_density_m, grad3_mom1_m, grad3_mom2_m, grad3_mom3_m, grad3_energy_m,  &
            normal_momentum, normal_grad, r, unorm_1, unorm_2, unorm_3

        real(rk),   allocatable, dimension(:,:) :: grid_velocity


        ! Interpolate interior solution to quadrature nodes
        density_m = worker%get_field('Density'   , 'value', 'face interior')
        mom1_m    = worker%get_field('Momentum-1', 'value', 'face interior')
        mom2_m    = worker%get_field('Momentum-2', 'value', 'face interior')
        mom3_m    = worker%get_field('Momentum-3', 'value', 'face interior')
        energy_m  = worker%get_field('Energy'    , 'value', 'face interior')


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


        ! Account for cylindrical. Get tangential momentum from angular momentum.
        r = worker%coordinate('1','boundary')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / r
        end if


        ! Initialize arrays
        density_bc = density_m
        mom1_bc    = mom1_m
        mom2_bc    = mom2_m
        mom3_bc    = mom3_m
        energy_bc  = energy_m


        ! Get unit normal vector
        unorm_1 = worker%unit_normal_ale(1)
        unorm_2 = worker%unit_normal_ale(2)
        unorm_3 = worker%unit_normal_ale(3)


        !   ANALYSIS
        !
        !   We want to reflect the normal part of the relative velocity.
        !   
        !   Relationship between absolute(U), relative(W), and ALE velocity(VG)
        !       U = W + VG
        !
        !   The condition we would like to impose is
        !       W_bc = W_m - 2*(W_m dot n)*n
        !
        !   Substituting in W = U - VG
        !       U_bc - VG = U_m - VG - 2*[(U_m - VG) dot n]*n
        !
        !   Simplifying gives
        !       U_bc = U_m - 2*[(U_m - VG) dot n]*n
        !


        ! Dot relative momentum with normal vector
        grid_velocity = worker%get_grid_velocity_face('face interior')
        normal_momentum = (mom1_m-grid_velocity(:,1)*density_m)*unorm_1 + &
                          (mom2_m-grid_velocity(:,2)*density_m)*unorm_2 + &
                          (mom3_m-grid_velocity(:,3)*density_m)*unorm_3

        ! Subtract relative normal momentum from relative momentum
        mom1_bc = mom1_m  -  normal_momentum*unorm_1
        mom2_bc = mom2_m  -  normal_momentum*unorm_2
        mom3_bc = mom3_m  -  normal_momentum*unorm_3


        ! Account for cylindrical. Convert tangential momentum back to angular momentum.
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_bc = mom2_bc * r
        end if


        ! Store boundary condition state
        call worker%store_bc_state('Density'   , density_bc, 'value')
        call worker%store_bc_state('Momentum-1', mom1_bc,    'value')
        call worker%store_bc_state('Momentum-2', mom2_bc,    'value')
        call worker%store_bc_state('Momentum-3', mom3_bc,    'value')
        call worker%store_bc_state('Energy'    , energy_bc,  'value')


        
        !call worker%store_bc_state('Density'   , ZERO*grad1_density_m, 'grad1')
        !call worker%store_bc_state('Density'   , ZERO*grad2_density_m, 'grad2')
        !call worker%store_bc_state('Density'   , ZERO*grad3_density_m, 'grad3')
        !                                        
        !call worker%store_bc_state('Momentum-1', ZERO*grad1_mom1_m,    'grad1')
        !call worker%store_bc_state('Momentum-1', ZERO*grad2_mom1_m,    'grad2')
        !call worker%store_bc_state('Momentum-1', ZERO*grad3_mom1_m,    'grad3')
        !                                        
        !call worker%store_bc_state('Momentum-2', ZERO*grad1_mom2_m,    'grad1')
        !call worker%store_bc_state('Momentum-2', ZERO*grad2_mom2_m,    'grad2')
        !call worker%store_bc_state('Momentum-2', ZERO*grad3_mom2_m,    'grad3')
        !                                        
        !call worker%store_bc_state('Momentum-3', ZERO*grad1_mom3_m,    'grad1')
        !call worker%store_bc_state('Momentum-3', ZERO*grad2_mom3_m,    'grad2')
        !call worker%store_bc_state('Momentum-3', ZERO*grad3_mom3_m,    'grad3')
        !                                        
        !call worker%store_bc_state('Energy'    , ZERO*grad1_energy_m,  'grad1')
        !call worker%store_bc_state('Energy'    , ZERO*grad2_energy_m,  'grad2')
        !call worker%store_bc_state('Energy'    , ZERO*grad3_energy_m,  'grad3')
        
        normal_grad = grad1_density_m*unorm_1 + grad2_density_m*unorm_2 + grad3_density_m*unorm_3
        grad1_density_m = grad1_density_m-normal_grad*unorm_1
        grad2_density_m = grad2_density_m-normal_grad*unorm_2
        grad3_density_m = grad3_density_m-normal_grad*unorm_3
        call worker%store_bc_state('Density'   , grad1_density_m, 'grad1')
        call worker%store_bc_state('Density'   , grad2_density_m, 'grad2')
        call worker%store_bc_state('Density'   , grad3_density_m, 'grad3')
        
        normal_grad = grad1_mom1_m*unorm_1 + grad2_mom1_m*unorm_2 + grad3_mom1_m*unorm_3
        grad1_mom1_m = grad1_mom1_m-normal_grad*unorm_1
        grad2_mom1_m = grad2_mom1_m-normal_grad*unorm_2
        grad3_mom1_m = grad3_mom1_m-normal_grad*unorm_3
        call worker%store_bc_state('Momentum-1', grad1_mom1_m,    'grad1')
        call worker%store_bc_state('Momentum-1', grad2_mom1_m,    'grad2')
        call worker%store_bc_state('Momentum-1', grad3_mom1_m,    'grad3')
                                                
        normal_grad = grad1_mom2_m*unorm_1 + grad2_mom2_m*unorm_2 + grad3_mom2_m*unorm_3
        grad1_mom2_m = grad1_mom2_m-normal_grad*unorm_1
        grad2_mom2_m = grad2_mom2_m-normal_grad*unorm_2
        grad3_mom2_m = grad3_mom2_m-normal_grad*unorm_3
        call worker%store_bc_state('Momentum-2', grad1_mom2_m,    'grad1')
        call worker%store_bc_state('Momentum-2', grad2_mom2_m,    'grad2')
        call worker%store_bc_state('Momentum-2', grad3_mom2_m,    'grad3')
                                                
        normal_grad = grad1_mom3_m*unorm_1 + grad2_mom3_m*unorm_2 + grad3_mom3_m*unorm_3
        grad1_mom3_m = grad1_mom3_m-normal_grad*unorm_1
        grad2_mom3_m = grad2_mom3_m-normal_grad*unorm_2
        grad3_mom3_m = grad3_mom3_m-normal_grad*unorm_3
        call worker%store_bc_state('Momentum-3', grad1_mom3_m,    'grad1')
        call worker%store_bc_state('Momentum-3', grad2_mom3_m,    'grad2')
        call worker%store_bc_state('Momentum-3', grad3_mom3_m,    'grad3')
                                                
        normal_grad = grad1_energy_m*unorm_1 + grad2_energy_m*unorm_2 + grad3_energy_m*unorm_3
        grad1_energy_m = grad1_energy_m-normal_grad*unorm_1
        grad2_energy_m = grad2_energy_m-normal_grad*unorm_2
        grad3_energy_m = grad3_energy_m-normal_grad*unorm_3
        call worker%store_bc_state('Energy'    , grad1_energy_m,  'grad1')
        call worker%store_bc_state('Energy'    , grad2_energy_m,  'grad2')
        call worker%store_bc_state('Energy'    , grad3_energy_m,  'grad3')


    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_symmetry
