module HP_bc
    use mod_kinds,          only: ik, rk
    use mod_constants,      only: ZERO, ONE, TWO
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
    type, public, extends(operator_t) :: HP_bc_t


    contains

        procedure   :: init
        procedure   :: compute

    end type HP_bc_t
    !*******************************************************************************************

    real(rk) :: p_param = 4._rk



contains






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(HP_bc_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Hyperbolized Poisson BC Operator')

        !
        ! Set operator type
        !
        !call self%set_operator_type('BC Diffusive Operator')
        call self%set_operator_type('BC Advective Operator')

        !
        ! Set operator equations
        !
        call self%add_primary_field('u')
        call self%add_primary_field('p')
        call self%add_primary_field('q')
        call self%add_primary_field('r')

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
        class(HP_bc_t),    intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            u, p, q, r, &
            flux_1, flux_2, flux_3, sumsqr, mag


        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        u  = worker%get_field('u','value','boundary')
        p  = worker%get_field('p','value','boundary')
        q  = worker%get_field('q','value','boundary')
        r  = worker%get_field('r','value','boundary')




        ! Allocate derivatives
        flux_1 = u*ZERO
        flux_2 = u*ZERO
        flux_3 = u*ZERO

        
        !
        ! u-equation
        !
        sumsqr = p*p + q*q + r*r
        if (abs(p_param-2._rk) > 1.e-8_rk) then
            mag = sumsqr**((p_param-TWO)/TWO)
        else
            mag = sumsqr
            mag = ONE
        end if

        flux_1 = mag*p
        flux_2 = mag*q
        flux_3 = mag*r

        call worker%integrate_boundary_condition('u','Advection',flux_1,flux_2,flux_3)

        !
        ! p-equation
        !
        flux_1 = u
        flux_2 = ZERO
        flux_3 = ZERO

        call worker%integrate_boundary_condition('p','Advection',flux_1,flux_2,flux_3)


        !
        ! q-equation
        !
        flux_1 = ZERO
        flux_2 = u
        flux_3 = ZERO

        call worker%integrate_boundary_condition('q','Advection',flux_1,flux_2,flux_3)


        !
        ! r-equation
        !
        flux_1 = ZERO
        flux_2 = ZERO
        flux_3 = u

        call worker%integrate_boundary_condition('r','Advection',flux_1,flux_2,flux_3)

    end subroutine compute
    !**********************************************************************************************










end module HP_bc
