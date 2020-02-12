module fcn_pressure_pulse
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, EIGHT, PI
    use type_point,     only: point_t
    use type_point_ad,  only: point_ad_t
    use type_function,  only: function_t
    use DNAD_D
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
    type, extends(function_t), public :: pressure_pulse_f

        real(rk)    :: x0 = 0.0_rk
        real(rk)    :: y0 = 0.0_rk
        real(rk)    :: z0 = 0.0_rk

        real(rk)    :: uinf   = ZERO
        real(rk)    :: vinf   = ZERO
        real(rk)    :: winf   = ZERO
        real(rk)    :: rhoinf = ONE
        real(rk)    :: pinf   = ONE/1.4_rk

        integer(ik) :: ivar

    contains

        procedure   :: init
        procedure   :: compute

    end type pressure_pulse_f
    !********************************************************************************



contains




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(pressure_pulse_f),  intent(inout)  :: self

        ! Set function name
        call self%set_name("pressure_pulse")

        call self%add_option('ivar', 1._rk)


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
        class(pressure_pulse_f),    intent(inout)   :: self
        real(rk),                   intent(in)      :: time
        type(point_ad_t),           intent(in)      :: coord

        type(AD_D)  :: x, y, z, r, &
                       u, v, w, &
                       rho, p, drho, dp, val

        integer(ik) :: ivar
        real(rk)    :: b, gam

        x   = coord%c1_
        y   = coord%c2_
        z   = coord%c3_
        b   = 0.01_rk
        gam = 1.4_rk


        !r = x
        !r = sqrt(x*x + y*y)
        r = sqrt(x*x + z*z)
        drho = 0.1_rk * exp(-(log(TWO)/b)*r*r)
        dp   = drho/gam


        rho = self%rhoinf + drho
        p   = self%pinf   + dp
        u   = ZERO
        v   = ZERO
        w   = ZERO




        ivar = nint(self%get_option_value('ivar'))
        select case (ivar)
            ! RHO
            case (1)
                val = rho

            ! RHO-U
            case (2)
                val = ZERO

            ! RHO-V
            case (3)
                val = ZERO

            ! RHO-W
            case (4)
                val = ZERO

            ! RHO-E
            case (5)
                val = p/(gam-ONE)  +  HALF*rho*(u*u + v*v + w*w)

            ! Pressure_TEMP
            case (6)
                val = p

            ! Grad1 P
            case (7)
                val = -13.8629*x*exp(-69.3147_rk*r*r)

            ! Grad2 P
            case (8)
                val = ZERO

            ! Grad3 P
            case (9)
                val = -13.8629*z*exp(-69.3147_rk*r*r)

        end select

    end function compute
    !**********************************************************************************


end module fcn_pressure_pulse
