module SA_LaxFriedrichs_operator
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF, ME, NEIGHBOR

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------
    type, extends(operator_t), public :: SA_LaxFriedrichs_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type SA_LaxFriedrichs_operator_t
    !**************************************************************************

contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(SA_LaxFriedrichs_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name("Scalar Advection LaxFriedrichs Operator")

        ! Set operator type
        call self%set_operator_type("Boundary Advective Operator")

        ! Set operator equations
        call self%set_equation("u")

    end subroutine init
    !********************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(SA_LaxFriedrichs_operator_t), intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        integer(ik)              :: iu
        real(rk)                 :: cx, cy, cz

        type(AD_D), dimension(:), allocatable   ::  &
            u_m,  u_p,                              &
            cx_m, cy_m, cz_m,                       &
            cx_p, cy_p, cz_p,                       &
            flux_x, flux_y, flux_z, integrand

        real(rk),   dimension(:), allocatable   ::  &
            normx, normy, normz, unormx, unormy, unormz


        !
        ! Get integer data
        !
        iu = prop%get_equation_index("u")


        !
        ! Interpolate solution to quadrature nodes
        !
        u_m = worker%interpolate(iu, 'value', ME)
        u_p = worker%interpolate(iu, 'value', NEIGHBOR)


        !
        ! Get model coefficients
        !
        cx_m = prop%scalar%compute_cx(u_m)
        cy_m = prop%scalar%compute_cy(u_m)
        cz_m = prop%scalar%compute_cz(u_m)
        cx_p = prop%scalar%compute_cx(u_p)
        cy_p = prop%scalar%compute_cy(u_p)
        cz_p = prop%scalar%compute_cz(u_p)


        !
        ! Get normal vector
        !
        normx  = worker%normal(1)
        normy  = worker%normal(2)
        normz  = worker%normal(3)

        unormx = worker%unit_normal(1)
        unormy = worker%unit_normal(2)
        unormz = worker%unit_normal(3)


        !
        ! Compute boundary upwind flux
        !
        flux_x = HALF*(cx_m*u_m - cx_p*u_p)
        flux_y = HALF*(cy_m*u_m - cy_p*u_p)
        flux_z = HALF*(cz_m*u_m - cz_p*u_p)

        integrand = flux_x*normx*unormx + flux_y*normy*unormy + flux_z*normz*unormz



        !
        ! Integrate flux
        !
        call worker%integrate_boundary(iu,integrand)


    end subroutine compute
    !*******************************************************************************************************









end module SA_LaxFriedrichs_operator
