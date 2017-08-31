module MMD_bc_operator
    use mod_kinds,          only: ik, rk
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    use ieee_arithmetic
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    type, public, extends(operator_t) :: MMD_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type MMD_bc_operator_t
    !*******************************************************************************************




contains






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(MMD_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Mesh Motion Diffusion BC Operator')

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
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!  @param[in]      face    face_info_t containing indices on location and misc information
    !!  @param[in]      flux    function_into_t containing info on the function being computed
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(MMD_bc_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            grad1_u1, grad2_u1, grad3_u1,           &
            grad1_u2, grad2_u2, grad3_u2,           &
            grad1_u3, grad2_u3, grad3_u3,           &
            flux_1,  flux_2,  flux_3, mu


        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        grad1_u1 = worker%get_field('grid_displacement1','grad1', 'boundary')
        grad2_u1 = worker%get_field('grid_displacement1','grad2', 'boundary')
        grad3_u1 = worker%get_field('grid_displacement1','grad3', 'boundary')

        grad1_u2 = worker%get_field('grid_displacement2','grad1', 'boundary')
        grad2_u2 = worker%get_field('grid_displacement2','grad2', 'boundary')
        grad3_u2 = worker%get_field('grid_displacement2','grad3', 'boundary')

        grad1_u3 = worker%get_field('grid_displacement3','grad1', 'boundary')
        grad2_u3 = worker%get_field('grid_displacement3','grad2', 'boundary')
        grad3_u3 = worker%get_field('grid_displacement3','grad3', 'boundary')


        !
        ! Compute scalar coefficient
        !
        mu = worker%get_field('Mesh Motion Diffusion Coefficient', 'value', 'boundary')



        !=================================================
        ! GD1 flux
        !=================================================
        flux_1 = -mu*grad1_u1
        flux_2 = -mu*grad2_u1
        flux_3 = -mu*grad3_u1

        call worker%integrate_boundary_condition('grid_displacement1','Diffusion',flux_1,flux_2,flux_3)


        !=================================================
        ! GD2 flux
        !=================================================
        flux_1 = -mu*grad1_u2
        flux_2 = -mu*grad2_u2
        flux_3 = -mu*grad3_u2

        call worker%integrate_boundary_condition('grid_displacement2','Diffusion',flux_1,flux_2,flux_3)


        !=================================================
        ! GD3 flux
        !=================================================
        flux_1 = -mu*grad1_u3
        flux_2 = -mu*grad2_u3
        flux_3 = -mu*grad3_u3

        call worker%integrate_boundary_condition('grid_displacement3','Diffusion',flux_1,flux_2,flux_3)


    end subroutine compute
    !**********************************************************************************************










end module MMD_bc_operator
