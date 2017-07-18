module fcn_convecting_vortex
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, EIGHT, PI
    use type_point,     only: point_t
    use type_function,  only: function_t
    implicit none
    private








    !>  This function defines an analytical solution to the 2D Euler equations.
    !!  It may be used to provide initial conditions for an unsteady simulation,
    !!  as well as to compute numerical errors.
    !!
    !!  This solution is implemented as in Fidkowski 2016
    !!
    !!  @author Eric Wolf
    !!  @date   2/15/2017
    !!
    !!
    !-------------------------------------------------------------------------------
    type, extends(function_t), public :: convecting_vortex_f
        private


        real(rk)    :: rho = ONE
        real(rk)    :: p   = ONE


        real(rk)    :: gam  = 1.4_rk
        real(rk)    :: vortex_strength = 0.3_rk !vortex strength
        real(rk)    :: rc = 1.5_rk !vortex size

        real(rk)    :: x0 = 5.0_rk
        real(rk)    :: y0 = 5.0_rk
        real(rk)    :: z0 = 0.0_rk

        real(rk)    :: uinf = TWO/sqrt(FIVE) 
        real(rk)    :: vinf = ONE/sqrt(FIVE)
        real(rk)    :: winf = ZERO


        real(rk)    :: rhoinf = ONE
        real(rk)    :: pinf = 20.0_rk/7.0_rk



        integer(ik) :: ivar

    contains

        procedure   :: init
        procedure   :: compute

    end type convecting_vortex_f
    !********************************************************************************



contains




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(convecting_vortex_f),  intent(inout)  :: self

        !
        ! Set function name
        !
        !self%name = "isentropic vortex  ::  weeeee!"
        call self%set_name("convecting_vortex")

        call self%add_option('ivar', 1._rk)

        !
        ! Set function options to default settings
        !
!        call self%dict%set('uinf',2._rk/sqrt(5._rk))
!        call self%dict%set('vinf',1._rk/sqrt(5._rk))
!        call self%dict%set('winf',0._rk)
!        call self%dict%set('Minf',1._rk*sqrt(
!        call self%dict%set('beta',0._rk)
!        call self%dict%set('xo',1._rk)
!        call self%dict%set('yo',1._rk)
!        call self%dict%set('zo',1._rk)
!

    end subroutine init
    !*************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    impure elemental function compute(self,time,coord) result(val)
        class(convecting_vortex_f),     intent(inout)  :: self
        real(rk),                       intent(in)  :: time
        type(point_t),                  intent(in)  :: coord

        real(rk)                                    :: val
        integer(ik)                                 :: ivar

        real(rk)    :: x,   y,   z, &
                        dx, dy, &
                       du, dv, u, v, w, &
                       gam, Minf, r, speedinf, rho, p, &
                       f0, f1, f2

        x = coord%c1_
        y = coord%c2_
        z = coord%c3_

        gam = self%gam
        speedinf = sqrt(self%uinf**TWO + self%vinf**TWO)
        Minf = speedinf*sqrt(self%rhoinf/(gam*self%pinf))

        dx = x - self%x0 - self%uinf*time
        dy = y - self%y0 - self%vinf*time
        r = sqrt(dx**TWO + dy**TWO)

        f0 = ONE - r**TWO/self%rc**TWO
        f1 = ONE - self%vortex_strength**TWO*(gam-ONE)*Minf**TWO*exp(f0)/(8._rk*PI**TWO)
        f2 = speedinf*self%vortex_strength*exp(f0/TWO)/(TWO*PI*self%rc)

        u = self%uinf - f2*dy
        v = self%vinf + f2*dx
        w  = ZERO

        rho = self%rhoinf*f1**(ONE/(gam-ONE))
        p = self%pinf*f1**(gam/(gam-ONE))






        ivar = nint(self%get_option_value('ivar'))
        select case (ivar)
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


end module fcn_convecting_vortex
