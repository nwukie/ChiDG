!>  Mesh Motion with a Linear Elasticity Model
!!  Reference: Yang and Mavripilis, 2009
!!
!!
!!
!-----------------------------------------------------------------------------------------------------------------------


module MMLE_boundary_operator
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: ZERO,ONE,TWO,HALF
    use type_operator,              only: operator_t
    use type_chidg_worker,          only: chidg_worker_t
    use type_properties,            only: properties_t
    use DNAD_D
    use ieee_arithmetic
    implicit none

    private



    !>
    !!
    !!  @author Eric Wolf
    !!  @date 5/16/2017
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: MMLE_boundary_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type MMLE_boundary_operator_t
    !********************************************************************************



contains


    !>
    !!
    !!  @author Eric Wolf
    !!  @date   5/16/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(MMLE_boundary_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Mesh Motion Linear Elasticity Boundary Average Operator")

        !
        ! Set operator type
        !
        call self%set_operator_type("Boundary Diffusive Operator")

        !
        ! Set operator equations
        !
        call self%add_primary_field("grid_displacement1")
        call self%add_primary_field("grid_displacement2")
        call self%add_primary_field("grid_displacement3")

    end subroutine init
    !********************************************************************************







    !>  Compute the diffusive boundary flux for scalar linear diffusion.
    !!
    !!  @author Eric Wolf
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(MMLE_boundary_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop



        type(AD_D), allocatable, dimension(:)   ::  &
            grad1_u1_m, grad2_u1_m, grad3_u1_m,     &
            grad1_u1_p, grad2_u1_p, grad3_u1_p,     &
            grad1_u2_m, grad2_u2_m, grad3_u2_m,     &
            grad1_u2_p, grad2_u2_p, grad3_u2_p,     &
            grad1_u3_m, grad2_u3_m, grad3_u3_m,     &
            grad1_u3_p, grad2_u3_p, grad3_u3_p,     &
            strain_11_m, strain_22_m, strain_33_m,  &
            strain_12_m, strain_31_m, strain_23_m,  &
            stress_11_m, stress_22_m, stress_33_m,  &
            stress_12_m, stress_31_m, stress_23_m,  &
            strain_11_p, strain_22_p, strain_33_p,  &
            strain_12_p, strain_31_p, strain_23_p,  &
            stress_11_p, stress_22_p, stress_33_p,  &
            stress_12_p, stress_31_p, stress_23_p,  &
            elasticity_modulus_m, alpha_param_m,    &
            elasticity_modulus_p, alpha_param_p,    &
            poisson_ratio_m, poisson_ratio_p,       &
            mu_m, mu_p,                             &
            flux_1_m, flux_2_m, flux_3_m,           &
            flux_1_p, flux_2_p, flux_3_p
             

        !
        ! Interpolate solution to quadrature nodes
        !
        grad1_u1_m = worker%get_field('grid_displacement1', 'grad1', 'face interior')
        grad2_u1_m = worker%get_field('grid_displacement1', 'grad2', 'face interior')
        grad3_u1_m = worker%get_field('grid_displacement1', 'grad3', 'face interior')


        grad1_u1_p = worker%get_field('grid_displacement1', 'grad1', 'face exterior')
        grad2_u1_p = worker%get_field('grid_displacement1', 'grad2', 'face exterior')
        grad3_u1_p = worker%get_field('grid_displacement1', 'grad3', 'face exterior')

        grad1_u2_m = worker%get_field('grid_displacement2', 'grad1', 'face interior')
        grad2_u2_m = worker%get_field('grid_displacement2', 'grad2', 'face interior')
        grad3_u2_m = worker%get_field('grid_displacement2', 'grad3', 'face interior')


        grad1_u2_p = worker%get_field('grid_displacement2', 'grad1', 'face exterior')
        grad2_u2_p = worker%get_field('grid_displacement2', 'grad2', 'face exterior')
        grad3_u2_p = worker%get_field('grid_displacement2', 'grad3', 'face exterior')

        grad1_u3_m = worker%get_field('grid_displacement3', 'grad1', 'face interior')
        grad2_u3_m = worker%get_field('grid_displacement3', 'grad2', 'face interior')
        grad3_u3_m = worker%get_field('grid_displacement3', 'grad3', 'face interior')


        grad1_u3_p = worker%get_field('grid_displacement3', 'grad1', 'face exterior')
        grad2_u3_p = worker%get_field('grid_displacement3', 'grad2', 'face exterior')
        grad3_u3_p = worker%get_field('grid_displacement3', 'grad3', 'face exterior')


        strain_11_m = grad1_u1_m
        strain_22_m = grad2_u2_m
        strain_33_m = grad3_u3_m

        strain_12_m = (grad2_u1_m + grad1_u2_m)
        strain_23_m = (grad3_u2_m + grad2_u3_m)
        strain_31_m = (grad3_u1_m + grad1_u3_m)


        !
        ! Compute scalar coefficient
        ! 
        poisson_ratio_m      = worker%get_field('Mesh Motion Linear Elasticity Poisson Ratio', 'value', 'face interior')
        elasticity_modulus_m = worker%get_field('Mesh Motion Linear Elasticity Modulus',       'value', 'face interior')


        !Use negative Poisson ratio to maintain element aspect ratio
        !elasticity_modulus_m = ONE
        alpha_param_m = elasticity_modulus_m/((ONE+poisson_ratio_m)*(ONE-TWO*poisson_ratio_m))


        stress_11_m = alpha_param_m*(&
            (ONE-poisson_ratio_m)*strain_11_m + poisson_ratio_m*strain_22_m + poisson_ratio_m*strain_33_m)
        stress_22_m = alpha_param_m*(&
            poisson_ratio_m*strain_11_m + (ONE-poisson_ratio_m)*strain_22_m + poisson_ratio_m*strain_33_m)
        stress_33_m = alpha_param_m*(&
            poisson_ratio_m*strain_11_m + poisson_ratio_m*strain_22_m + (ONE-poisson_ratio_m)*strain_33_m)

        stress_12_m = alpha_param_m*(HALF-poisson_ratio_m)*strain_12_m
        stress_23_m = alpha_param_m*(HALF-poisson_ratio_m)*strain_23_m
        stress_31_m = alpha_param_m*(HALF-poisson_ratio_m)*strain_31_m
        


        strain_11_p = grad1_u1_p
        strain_22_p = grad2_u2_p
        strain_33_p = grad3_u3_p

        strain_12_p = (grad2_u1_p + grad1_u2_p)
        strain_23_p = (grad3_u2_p + grad2_u3_p)
        strain_31_p = (grad3_u1_p + grad1_u3_p)

        !!
        !! Compute scalar coefficient
        !! 

        poisson_ratio_p = worker%get_field('Mesh Motion Linear Elasticity Poisson Ratio', 'value','face exterior')
        elasticity_modulus_p  = worker%get_field('Mesh Motion Linear Elasticity Modulus', 'value','face exterior')

        !elasticity_modulus_p = ONE
        alpha_param_p = elasticity_modulus_p/((ONE+poisson_ratio_p)*(ONE-TWO*poisson_ratio_p))


        stress_11_p = alpha_param_p*(&
            (ONE-poisson_ratio_p)*strain_11_p + poisson_ratio_p*strain_22_p + poisson_ratio_p*strain_33_p)
        stress_22_p = alpha_param_p*(&
            poisson_ratio_p*strain_11_p + (ONE-poisson_ratio_p)*strain_22_p + poisson_ratio_p*strain_33_p)
        stress_33_p = alpha_param_p*(&
            poisson_ratio_p*strain_11_p + poisson_ratio_p*strain_22_p + (ONE-poisson_ratio_p)*strain_33_p)

        stress_12_p = alpha_param_p*(HALF-poisson_ratio_p)*strain_12_p
        stress_23_p = alpha_param_p*(HALF-poisson_ratio_p)*strain_23_p
        stress_31_p = alpha_param_p*(HALF-poisson_ratio_p)*strain_31_p
        






        !====================================
        !               GD1
        !====================================
        flux_1_m = -stress_11_m
        flux_2_m = -stress_12_m
        flux_3_m = -stress_31_m

        flux_1_p = -stress_11_p
        flux_2_p = -stress_12_p
        flux_3_p = -stress_31_p

        call worker%integrate_boundary_average('grid_displacement1','Diffusion',    &
                                                flux_1_m, flux_2_m, flux_3_m,       &
                                                flux_1_p, flux_2_p, flux_3_p)


        !====================================
        !               GD1
        !====================================
        flux_1_m = -stress_12_m
        flux_2_m = -stress_22_m
        flux_3_m = -stress_23_m

        flux_1_p = -stress_12_p
        flux_2_p = -stress_22_p
        flux_3_p = -stress_23_p

        call worker%integrate_boundary_average('grid_displacement2','Diffusion',    &
                                                flux_1_m, flux_2_m, flux_3_m,       &
                                                flux_1_p, flux_2_p, flux_3_p)


        !====================================
        !               GD1
        !====================================
        flux_1_m = -stress_31_m
        flux_2_m = -stress_23_m
        flux_3_m = -stress_33_m

        flux_1_p = -stress_31_p
        flux_2_p = -stress_23_p
        flux_3_p = -stress_33_p

        call worker%integrate_boundary_average('grid_displacement3','Diffusion',    &
                                                flux_1_m, flux_2_m, flux_3_m,       &
                                                flux_1_p, flux_2_p, flux_3_p)


    end subroutine compute
    !**************************************************************************************************




end module MMLE_boundary_operator
