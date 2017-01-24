module SA_LaxFriedrichs_operator
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF
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
    !-----------------------------------------------------------------------------------------------
    type, extends(operator_t), public :: SA_LaxFriedrichs_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type SA_LaxFriedrichs_operator_t
    !***********************************************************************************************

contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine init(self)
        class(SA_LaxFriedrichs_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name("Scalar Advection LaxFriedrichs Operator")

        ! Set operator type
        call self%set_operator_type("Boundary Advective Operator")

        ! Set operator equations
        call self%add_primary_field("u")

    end subroutine init
    !***********************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(SA_LaxFriedrichs_operator_t), intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        integer(ik)              :: iu

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
        iu = prop%get_primary_field_index("u")


        !
        ! Interpolate solution to quadrature nodes
        !
        u_m    = worker%get_primary_field_face('u', 'value', 'face interior')
        u_p    = worker%get_primary_field_face('u', 'value', 'face exterior')


        !
        ! Get model coefficients
        !
        cx_m = worker%get_model_field_face('Scalar X-Advection Velocity', 'value', 'face interior')
        cy_m = worker%get_model_field_face('Scalar Y-Advection Velocity', 'value', 'face interior')
        cz_m = worker%get_model_field_face('Scalar Z-Advection Velocity', 'value', 'face interior')
        cx_p = worker%get_model_field_face('Scalar X-Advection Velocity', 'value', 'face exterior')
        cy_p = worker%get_model_field_face('Scalar Y-Advection Velocity', 'value', 'face exterior')
        cz_p = worker%get_model_field_face('Scalar Z-Advection Velocity', 'value', 'face exterior')
        

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
        call worker%integrate_boundary('u',integrand)


    end subroutine compute
    !*************************************************************************************************









end module SA_LaxFriedrichs_operator
