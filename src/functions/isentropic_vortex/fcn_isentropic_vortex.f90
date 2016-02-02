module fcn_isentropic_vortex
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, EIGHT, PI
    use type_point,     only: point_t
    use type_function,  only: function_t
    implicit none
    private








    !>  @TODO NEEDS FIXED
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    type, extends(function_t), public :: isentropic_vortex_f
        private


        real(rk)    :: rho = ONE
        real(rk)    :: p   = ONE


        real(rk)    :: gam  = 1.4_rk
        real(rk)    :: beta = FIVE
        real(rk)    :: xo
        real(rk)    :: yo
        real(rk)    :: zo

        real(rk)    :: uinf = ONE
        real(rk)    :: vinf = ONE
        real(rk)    :: winf = ZERO

        integer(ik) :: ivar

    contains

        procedure   :: init
        procedure   :: compute

    end type isentropic_vortex_f
    !********************************************************************************



contains




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(isentropic_vortex_f),  intent(inout)  :: self

        !
        ! Set function name
        !
        self%name = "isentropic vortex  ::  weeeee!"


        !
        ! Set function options to default settings
        !
        call self%dict%set('uinf',1._rk)
        call self%dict%set('vinf',0._rk)
        call self%dict%set('winf',0._rk)
        call self%dict%set('beta',0._rk)
        call self%dict%set('xo',1._rk)
        call self%dict%set('yo',1._rk)
        call self%dict%set('zo',1._rk)


    end subroutine init
    !*************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    elemental function compute(self,time,coord) result(val)
        class(isentropic_vortex_f),     intent(in)  :: self
        real(rk),                       intent(in)  :: time
        type(point_t),                  intent(in)  :: coord

        real(rk)                                    :: val

        real(rk)    :: x,   y,   z, &
                       du, dv, u, v, w, &
                       gam, beta, r, T, rho, p

        x = coord%c1_
        y = coord%c2_
        z = coord%c3_

        gam = self%gam
        beta = self%beta

        r = sqrt((x - self%xo)**TWO + (y - self%yo)**TWO)
        T = ONE - ((gam - ONE)*(beta**TWO)/(EIGHT*gam*PI**TWO))*exp(ONE - r**TWO)
        rho = T**(ONE/(gam-ONE))
        p   = rho*T

        du = (beta/(TWO*PI))*exp((ONE-r**TWO)/TWO)*(-(y-self%yo))
        dv = (beta/(TWO*PI))*exp((ONE-r**TWO)/TWO)*((x-self%xo))
        u  = self%uinf + du
        v  = self%vinf + dv
        w  = ZERO





        select case (self%ivar)
            ! RHO
            case (1)
                val = rho

            ! RHO-U
            case (2)
                val = rho*u

            ! RHO-V
            case (3)
                val = rho*v

            ! RHO-W
            case (4)
                val = rho*w

            ! RHO-E
            case (5)
                val = p/(gam-ONE)  +  HALF*rho*(u*u + v*v + w*w)

        end select

    end function compute
    !**********************************************************************************


end module fcn_isentropic_vortex
