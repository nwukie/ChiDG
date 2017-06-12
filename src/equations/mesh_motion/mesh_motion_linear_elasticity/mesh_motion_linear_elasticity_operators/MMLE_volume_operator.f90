!>  Mesh Motion with a Linear Elasticity Model
!!  Reference: Yang and Mavripilis, 2009
!!
!!
!!
!-----------------------------------------------------------------------------------------------------------------------

module MMLE_volume_operator
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    use ieee_arithmetic
    implicit none
    private

    !>
    !!
    !!  @author  Eric Wolf
    !!  @date   5/16/2017
    !!
    !!
    !!
    !-------------------------------------------------------------------------
    type, extends(operator_t), public :: MMLE_volume_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type MMLE_volume_operator_t
    !*************************************************************************

contains


    !>
    !!
    !!  @author Eric Wolf
    !!  @date   5/16/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(MMLE_volume_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name('Mesh Motion Linear Elasticity Volume Operator')

        ! Set operator type
        call self%set_operator_type('Volume Diffusive Operator')

        ! Set operator equations
        call self%add_primary_field('grid_displacement1')
        call self%add_primary_field('grid_displacement2')
        call self%add_primary_field('grid_displacement3')

    end subroutine init
    !********************************************************************************





    !>
    !!
    !!  @author Eric Wolf 
    !!  @date   5/16/2017
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(MMLE_volume_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:)   ::  &
            flux_1, flux_2, flux_3,  mu, &
            grad1_u1, grad2_u1, grad3_u1, &
            grad1_u2, grad2_u2, grad3_u2, &
            grad1_u3, grad2_u3, grad3_u3, &
            strain_11, strain_22, strain_33, &
            strain_12, strain_31, strain_23,&
            stress_11, stress_22, stress_33, &
            stress_12, stress_31, stress_23, &
            elasticity_modulus, alpha_param, poisson_ratio

        

        !
        ! Interpolate solution to quadrature nodes
        !
        grad1_u1 = worker%get_primary_field_element('grid_displacement1','grad1 + lift')
        grad2_u1 = worker%get_primary_field_element('grid_displacement1','grad2 + lift')
        grad3_u1 = worker%get_primary_field_element('grid_displacement1','grad3 + lift')

        grad1_u2 = worker%get_primary_field_element('grid_displacement2','grad1 + lift')
        grad2_u2 = worker%get_primary_field_element('grid_displacement2','grad2 + lift')
        grad3_u2 = worker%get_primary_field_element('grid_displacement2','grad3 + lift')

        grad1_u3 = worker%get_primary_field_element('grid_displacement3','grad1 + lift')
        grad2_u3 = worker%get_primary_field_element('grid_displacement3','grad2 + lift')
        grad3_u3 = worker%get_primary_field_element('grid_displacement3','grad3 + lift')



        strain_11 = grad1_u1
        strain_22 = grad2_u2
        strain_33 = grad3_u3

        strain_12 = (grad2_u1 + grad1_u2)
        strain_23 = (grad3_u2 + grad2_u3)
        strain_31 = (grad3_u1 + grad1_u3)

        !!
        !! Compute scalar coefficient
        !! 
        poisson_ratio = worker%get_model_field_element('Mesh Motion Linear Elasticity Poisson Ratio', 'value')
        elasticity_modulus = worker%get_model_field_element('Mesh Motion Linear Elasticity Modulus', 'value')


        !Use negative Poisson ratio to maintain element aspect ratio
        alpha_param = elasticity_modulus/((ONE+poisson_ratio)*(ONE-TWO*poisson_ratio))


        stress_11 = alpha_param*(&
            (ONE-poisson_ratio)*strain_11 + poisson_ratio*strain_22 + poisson_ratio*strain_33)
        stress_22 = alpha_param*(&
            poisson_ratio*strain_11 + (ONE-poisson_ratio)*strain_22 + poisson_ratio*strain_33)
        stress_33 = alpha_param*(&
            poisson_ratio*strain_11 + poisson_ratio*strain_22 + (ONE-poisson_ratio)*strain_33)

        stress_12 = alpha_param*(HALF-poisson_ratio)*strain_12
        stress_23 = alpha_param*(HALF-poisson_ratio)*strain_23
        stress_31 = alpha_param*(HALF-poisson_ratio)*strain_31
        


        
        !
        ! Compute volume flux at quadrature nodes
        !

        !GD1
        flux_1 = -stress_11
        flux_2 = -stress_12
        flux_3 = -stress_31


        !
        ! Integrate volume flux
        !
        call worker%integrate_volume('grid_displacement1',flux_1,flux_2,flux_3)

        !GD2
        flux_1 = -stress_12
        flux_2 = -stress_22
        flux_3 = -stress_23



        !
        ! Integrate volume flux
        !
        call worker%integrate_volume('grid_displacement2',flux_1,flux_2,flux_3)

        !GD3
        flux_1 = -stress_31
        flux_2 = -stress_23
        flux_3 = -stress_33



        !
        ! Integrate volume flux
        !
        call worker%integrate_volume('grid_displacement3',flux_1,flux_2,flux_3)





    end subroutine compute
    !****************************************************************************************************






end module MMLE_volume_operator
