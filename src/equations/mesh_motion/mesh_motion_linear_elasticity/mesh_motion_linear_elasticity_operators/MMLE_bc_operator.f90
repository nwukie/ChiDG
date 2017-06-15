!>  Mesh Motion with a Linear Elasticity Model
!!  Reference: Yang and Mavripilis, 2009
!!
!!
!!
!-----------------------------------------------------------------------------------------------------------------------


module MMLE_bc_operator
    use mod_kinds,          only: ik, rk
    use mod_constants,      only: ZERO, ONE, TWO, HALF
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    use ieee_arithmetic
    implicit none




    !>
    !!
    !!  @author Eric Wolf
    !!  @date   5/16/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    type, public, extends(operator_t) :: MMLE_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type MMLE_bc_operator_t
    !*******************************************************************************************




contains






    !>
    !!
    !!  @author Eric Wolf
    !!  @date   5/16/2017
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(MMLE_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Mesh Motion Linear Elasticity BC Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('BC Diffusive Operator')

        !
        ! Set operator equations
        !
        call self%add_primary_field('grid_displacement1')
        call self%add_primary_field('grid_displacement2')
        call self%add_primary_field('grid_displacement3')

    end subroutine init
    !********************************************************************************





    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Eric Wolf
    !!  @date   5/16/2017
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!  @param[in]      face    face_info_t containing indices on location and misc information
    !!  @param[in]      flux    function_into_t containing info on the function being computed
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(MMLE_bc_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            grad1_u1, grad2_u1, grad3_u1,              &
            grad1_u2, grad2_u2, grad3_u2,              &
            grad1_u3, grad2_u3, grad3_u3,              &
            flux_1,  flux_2,  flux_3,               &
            integrand, mu, &
            strain_11, strain_22, strain_33, &
            strain_12, strain_31, strain_23,&
            stress_11, stress_22, stress_33, &
            stress_12, stress_31, stress_23, &
            elasticity_modulus, alpha_param, poisson_ratio

        

        real(rk),   allocatable, dimension(:)   ::  &
            norm_1, norm_2, norm_3


        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        grad1_u1 = worker%get_primary_field_face('grid_displacement1','grad1 + lift', 'boundary')
        grad2_u1 = worker%get_primary_field_face('grid_displacement1','grad2 + lift', 'boundary')
        grad3_u1 = worker%get_primary_field_face('grid_displacement1','grad3 + lift', 'boundary')

        grad1_u2 = worker%get_primary_field_face('grid_displacement2','grad1 + lift', 'boundary')
        grad2_u2 = worker%get_primary_field_face('grid_displacement2','grad2 + lift', 'boundary')
        grad3_u2 = worker%get_primary_field_face('grid_displacement2','grad3 + lift', 'boundary')

        grad1_u3 = worker%get_primary_field_face('grid_displacement3','grad1 + lift', 'boundary')
        grad2_u3 = worker%get_primary_field_face('grid_displacement3','grad2 + lift', 'boundary')
        grad3_u3 = worker%get_primary_field_face('grid_displacement3','grad3 + lift', 'boundary')




        strain_11 = grad1_u1
        strain_22 = grad2_u2
        strain_33 = grad3_u3

        strain_12 = (grad2_u1 + grad1_u2)
        strain_23 = (grad3_u2 + grad2_u3)
        strain_31 = (grad3_u1 + grad1_u3)

        !!
        !! Compute scalar coefficient
        !! 
        poisson_ratio = worker%get_model_field_face('Mesh Motion Linear Elasticity Poisson Ratio', 'value', 'boundary')
        elasticity_modulus = worker%get_model_field_face('Mesh Motion Linear Elasticity Modulus', 'value', 'boundary')


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
 



        norm_1 = worker%normal(1)
        norm_2 = worker%normal(2)
        norm_3 = worker%normal(3)


        

        !=================================================
        ! GD1 flux
        !=================================================
        flux_1 = -stress_11
        flux_2 = -stress_12
        flux_3 = -stress_31


        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3


        if (any(ieee_is_nan(integrand(:)%x_ad_))) then
            print*, 'BC OP: ', integrand(:)%x_ad_
        end if


        call worker%integrate_boundary('grid_displacement1',integrand)


        !=================================================
        ! GD2 flux
        !=================================================

        flux_1 = -stress_12
        flux_2 = -stress_22
        flux_3 = -stress_23
        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3


        if (any(ieee_is_nan(integrand(:)%x_ad_))) then
            print*, 'BC OP: ', integrand(:)%x_ad_
        end if


        call worker%integrate_boundary('grid_displacement2',integrand)


        !=================================================
        ! GD3 flux
        !=================================================

        flux_1 = -stress_31
        flux_2 = -stress_23
        flux_3 = -stress_33
        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3


        if (any(ieee_is_nan(integrand(:)%x_ad_))) then
            print*, 'BC OP: ', integrand(:)%x_ad_
        end if


        call worker%integrate_boundary('grid_displacement3',integrand)


    end subroutine compute
    !**********************************************************************************************










end module MMLE_bc_operator
