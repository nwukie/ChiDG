module pgradtest_bc_operator
    use mod_kinds,          only: ik, rk
    use mod_constants,      only: ZERO, ONE
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
    type, public, extends(operator_t) :: pgradtest_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type pgradtest_bc_operator_t
    !*******************************************************************************************




contains






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(pgradtest_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('pgradtest BC Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('BC Diffusive Operator')

        !
        ! Set operator equations
        !
        call self%add_primary_field('Pressure')

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
        class(pgradtest_bc_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            grad1_pbc, grad2_pbc, grad3_pbc,        &
            flux_1,  flux_2,  flux_3, mu


        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        grad1_pbc = worker%get_field('Pressure', 'grad1', 'boundary')
        grad2_pbc = worker%get_field('Pressure', 'grad2', 'boundary')
        grad3_pbc = worker%get_field('Pressure', 'grad3', 'boundary')

        mu = grad1_pbc
        mu = ONE

        !=================================================
        ! Mass flux
        !=================================================
        flux_1 = mu*grad1_pbc
        flux_2 = mu*grad2_pbc
        flux_3 = mu*grad3_pbc

        call worker%integrate_boundary_condition('Pressure','Diffusion',flux_1,flux_2,flux_3)

    end subroutine compute
    !**********************************************************************************************










end module pgradtest_bc_operator
