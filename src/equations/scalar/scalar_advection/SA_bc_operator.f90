module SA_bc_operator
    use mod_kinds,          only: ik, rk
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    type, public, extends(operator_t) :: SA_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type SA_bc_operator_t
    !*******************************************************************************************




contains






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(SA_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Scalar Advection BC Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('BC Diffusive Operator')

        !
        ! Set operator equations
        !
        call self%add_primary_field('u')

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
        class(SA_bc_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            u, c1, c2, c3,                          &
            flux_1, flux_2, flux_3

        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        u  = worker%get_field('u', 'value', 'boundary')


        !
        ! Get model coefficients
        !
        c1 = worker%get_field('Scalar Advection Velocity-1', 'value', 'boundary')
        c2 = worker%get_field('Scalar Advection Velocity-2', 'value', 'boundary')
        c3 = worker%get_field('Scalar Advection Velocity-3', 'value', 'boundary')

        
        !
        ! Compute flux and integrate
        !
        flux_1 = c1*u
        flux_2 = c2*u
        flux_3 = c3*u

        call worker%integrate_boundary_condition('u','Advection',flux_1,flux_2,flux_3)


    end subroutine compute
    !**********************************************************************************************










end module SA_bc_operator
