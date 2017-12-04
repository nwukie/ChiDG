module DLA_bc_operator
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
    type, public, extends(operator_t) :: DLA_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type DLA_bc_operator_t
    !*******************************************************************************************




contains






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(DLA_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('DLA Advection BC Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('BC Advective Flux')

        !
        ! Set operator equations
        !
        call self%add_primary_field('u_a')
        call self%add_primary_field('u_b')

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
        class(DLA_bc_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            u_a, u_b, flux_1, flux_2, flux_3

        real(rk)    :: c1, c2, c3

        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        u_a = worker%get_field('u_a', 'value', 'boundary')
        u_b = worker%get_field('u_b', 'value', 'boundary')


        !
        ! Get model coefficients
        !
        c1 = 1._rk
        c2 = 0._rk
        c3 = 0._rk

        
        !
        ! Compute flux and integrate
        !
        flux_1 = c1*u_a
        flux_2 = c2*u_a
        flux_3 = c3*u_a

        call worker%integrate_boundary_condition('u_a','Advection',flux_1,flux_2,flux_3)


        !
        ! Compute flux and integrate
        !
        flux_1 = c1*u_b
        flux_2 = c2*u_b
        flux_3 = c3*u_b

        call worker%integrate_boundary_condition('u_b','Advection',flux_1,flux_2,flux_3)



    end subroutine compute
    !**********************************************************************************************










end module DLA_bc_operator
