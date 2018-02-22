module graddemo_P_bc_operator
    use mod_kinds,          only: ik, rk
    use mod_constants,      only: ZERO, ONE, HALF
    use mod_fluid,          only: gam, Rgas
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
    !---------------------------------------------------------------------------------------
    type, public, extends(operator_t) :: graddemo_P_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type graddemo_P_bc_operator_t
    !***************************************************************************************




contains






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine init(self)
        class(graddemo_P_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Graddemo P BC Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('BC Diffusive Operator')

        !
        ! Set operator equations
        !
        call self%add_primary_field('Pressure_TEMP')

    end subroutine init
    !**************************************************************************************





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
    !---------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(graddemo_P_bc_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            grad1_pbc, grad2_pbc, grad3_pbc,        &
            grad1_p,   grad2_p,   grad3_p,          &
            flux_1,    flux_2,    flux_3




        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        grad1_pbc = worker%get_field('Pressure_TEMP', 'grad1', 'boundary')
        grad2_pbc = worker%get_field('Pressure_TEMP', 'grad2', 'boundary')
        grad3_pbc = worker%get_field('Pressure_TEMP', 'grad3', 'boundary')

        grad1_p   = worker%get_field('Pressure Gradient - 1', 'value', 'boundary')
        grad2_p   = worker%get_field('Pressure Gradient - 2', 'value', 'boundary')
        grad3_p   = worker%get_field('Pressure Gradient - 3', 'value', 'boundary')


        flux_1 = grad1_pbc - grad1_p
        flux_2 = grad2_pbc - grad2_p
        flux_3 = grad3_pbc - grad3_p


        call worker%integrate_boundary_condition('Pressure_TEMP','Diffusion',flux_1,flux_2,flux_3)



    end subroutine compute
    !***************************************************************************************










end module graddemo_P_bc_operator
