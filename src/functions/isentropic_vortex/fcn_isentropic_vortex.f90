module fcn_isentropic_vortex
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
    !-------------------------------------------------------------------------------
    type, extends(function_t), public :: isentropic_vortex_f

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

        ! Set function name
        call self%set_name("isentropic_vortex")

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
        class(isentropic_vortex_f),     intent(inout)  :: self
        real(rk),                       intent(in)  :: time
        type(point_ad_t),               intent(in)  :: coord


        type(AD_D)  :: x, y, z, r, &
                       u, v, w, &
                       rho, p, drho, dp, theta, du, dv, val

        integer(ik) :: ivar
        real(rk)    :: gam, uinf, vinf, dumax, b

        call chidg_signal(FATAL,"isentropic_vortex: needs reimplemented for differentiated coordinates.")

!        x     = coord%c1_
!        y     = coord%c2_
!        z     = coord%c3_
!        r     = sqrt(x*x + y*y)
!!        theta = atan2(y,x)
!
!
!
!        !
!        ! Base parameters
!        !
!        gam = 1.4_rk
!        uinf = 0.5_rk
!        vinf = 0._rk
!        dumax = 0.5_rk*uinf
!        b = 0.2_rk
!
!        
!        !
!        ! Perturbations
!        !
!        drho = (ONE - HALF*(gam - ONE)*dumax*dumax*exp(ONE - (r*r/(b*b))))**(ONE/(gam-ONE))
!        dp   = (ONE/gam)*drho**(gam)
!        du   = -(dumax/b)*r*exp(HALF*(ONE - (r*r/(b*b)))) * sin(theta)
!        dv   =  (dumax/b)*r*exp(HALF*(ONE - (r*r/(b*b)))) * cos(theta)
!
!
!        !
!        ! Total quantities
!        !
!        rho = drho
!        p   = dp
!        u   = uinf + du
!        v   = vinf + dv
!        w   = ZERO
!
!
!
!
!        ivar = nint(self%get_option_value('ivar'))
!        select case (ivar)
!            ! RHO
!            case (1)
!                val = rho
!
!            ! RHO-U
!            case (2)
!                val = rho * u
!
!            ! RHO-V
!            case (3)
!                val = rho * v
!
!            ! RHO-W
!            case (4)
!                val = rho * w
!
!            ! RHO-E
!            case (5)
!                val = p/(gam-ONE)  +  HALF*rho*(u*u + v*v + w*w)
!
!        end select

    end function compute
    !**********************************************************************************


end module fcn_isentropic_vortex
