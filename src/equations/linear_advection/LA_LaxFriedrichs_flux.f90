module LA_LaxFriedrichs_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF, ME, NEIGHBOR

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D

!    use LA_properties,          only: LA_properties_t
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------
    type, extends(operator_t), public :: LA_LaxFriedrichs_flux_t


    contains

        procedure   :: init
        procedure   :: compute

    end type LA_LaxFriedrichs_flux_t
    !**************************************************************************

contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(LA_LaxFriedrichs_flux_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name("LA LaxFriedrichs Flux")

        ! Set operator type
        call self%set_operator_type("Boundary Advective Flux")

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
        class(LA_LaxFriedrichs_flux_t),     intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        integer(ik)              :: iu
        real(rk)                 :: cx, cy, cz

        type(AD_D), dimension(:), allocatable   ::  &
            u_l, u_r, flux_x, flux_y, flux_z, integrand

        real(rk),   dimension(:), allocatable   ::  &
            normx, normy, normz, unormx, unormy, unormz


        !
        ! Get integer data
        !
        iu = prop%get_equation_index("u")



        !
        ! Get equation set property data
        !
!        select type(prop)
!            type is (LA_properties_t)
!                cx = prop%c(1)
!                cy = prop%c(2)
!                cz = prop%c(3)
!        end select
        cx = 1._rk
        cy = 0._rk
        cz = 0._rk



        !
        ! Interpolate solution to quadrature nodes
        !
        u_r = worker%interpolate(iu, 'value', ME)
        u_l = worker%interpolate(iu, 'value', NEIGHBOR)


        normx  = worker%normal(1)
        normy  = worker%normal(2)
        normz  = worker%normal(3)

        unormx = worker%unit_normal(1)
        unormy = worker%unit_normal(2)
        unormz = worker%unit_normal(3)


        !
        ! Compute boundary upwind flux. [  Significant: (u_r - u_l) not (u_l - u_r)  ]
        !
        flux_x = (cx * (u_r - u_l)/TWO ) * normx * unormx
        flux_y = (cy * (u_r - u_l)/TWO ) * normy * unormy
        flux_z = (cz * (u_r - u_l)/TWO ) * normz * unormz

        integrand = flux_x + flux_y + flux_z



        !
        ! Integrate flux
        !
        call worker%integrate_boundary(iu,integrand)


    end subroutine compute
    !*******************************************************************************************************









end module LA_LaxFriedrichs_flux
