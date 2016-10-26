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
            dudx_m,dudy_m,dudz_m,                   &
            dudx_p,dudy_p,dudz_p,                   &
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
        u_m    = worker%get_face_variable(iu, 'value', ME)
        dudx_m = worker%get_face_variable(iu, 'ddx'  , ME)
        dudy_m = worker%get_face_variable(iu, 'ddy'  , ME)
        dudz_m = worker%get_face_variable(iu, 'ddz'  , ME)

        u_p    = worker%get_face_variable(iu, 'value', NEIGHBOR)
        dudx_p = worker%get_face_variable(iu, 'ddx'  , NEIGHBOR)
        dudy_p = worker%get_face_variable(iu, 'ddy'  , NEIGHBOR)
        dudz_p = worker%get_face_variable(iu, 'ddz'  , NEIGHBOR)


        !
        ! Get model coefficients
        !
        cx_m = prop%scalar%compute_cx(u_m,dudx_m,dudy_m,dudz_m)
        cy_m = prop%scalar%compute_cy(u_m,dudx_m,dudy_m,dudz_m)
        cz_m = prop%scalar%compute_cz(u_m,dudx_m,dudy_m,dudz_m)
        cx_p = prop%scalar%compute_cx(u_p,dudx_p,dudy_p,dudz_p)
        cy_p = prop%scalar%compute_cy(u_p,dudx_p,dudy_p,dudz_p)
        cz_p = prop%scalar%compute_cz(u_p,dudx_p,dudy_p,dudz_p)


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
        flux_x = max(abs(cx_m),abs(cx_p))*HALF*(u_m - u_p)
        flux_y = max(abs(cy_m),abs(cy_p))*HALF*(u_m - u_p)
        flux_z = max(abs(cz_m),abs(cz_p))*HALF*(u_m - u_p)

        integrand = flux_x*normx*unormx + flux_y*normy*unormy + flux_z*normz*unormz



        !
        ! Integrate flux
        !
        call worker%integrate_boundary(iu,integrand)


    end subroutine compute
    !*******************************************************************************************************









end module SA_LaxFriedrichs_operator
