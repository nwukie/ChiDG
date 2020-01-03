module fcn_cantilevered_beam
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FIVE, SIX, EIGHT, PI
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
    type, extends(function_t), public :: cantilevered_beam_f
        private

    contains

        procedure   :: init
        procedure   :: compute

    end type cantilevered_beam_f
    !********************************************************************************



contains




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(cantilevered_beam_f),  intent(inout)  :: self

        !
        ! Set function name
        !
        !self%name = "isentropic vortex  ::  weeeee!"
        call self%set_name("cantilevered_beam")

        call self%add_option('ivar', 1._rk)
        call self%add_option('Poisson Ratio', 0.0_rk)
        call self%add_option('Elasticity Modulus', 1._rk)
        call self%add_option('P',1._rk)
        call self%add_option('L',1._rk)
        call self%add_option('h',1._rk)
        call self%add_option('w',1._rk)

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
        class(cantilevered_beam_f),     intent(inout)  :: self
        real(rk),                       intent(in)  :: time
        type(point_t),                  intent(in)  :: coord

        real(rk)                                    :: val
        integer(ik)                                 :: ivar

        real(rk)    :: x,   y,   z, &
                        P, L, h, w, &
                    poisson_ratio, elasticity_modulus, &
                    alpha, shear_modulus, gd1, gd2, gd3


        P = self%get_option_value('P')
        L = self%get_option_value('L')
        h = self%get_option_value('h')
        w = self%get_option_value('w')

        x = coord%c1_ 
        y = coord%c2_ - h/TWO
        z = coord%c3_


        poisson_ratio = self%get_option_value('Poisson Ratio')
        elasticity_modulus = self%get_option_value('Elasticity Modulus')


        alpha = ONE/(ONE+poisson_ratio)
        shear_modulus = elasticity_modulus/(TWO*(ONE+poisson_ratio))


        gd1 = HALF*(TWO*P*y**THREE/(h**THREE*w)-THREE*HALF*y*P/(h*w)+ &
            alpha*(-SIX*x**TWO*P*y/(h**THREE*w) + &
            SIX*TWO*L*P*y*x/(h**THREE*w)+ &
            TWO*P*y**THREE/(h**THREE*w) - &
            HALF*(-TWO*P+alpha*P)*y/(w*alpha*h)))/shear_modulus
        
        gd2 = HALF*(SIX*x*P*y**TWO/(h**THREE*w)-SIX*x*P*L*y**TWO/(h**THREE*w)-&
            THREE*HALF*x*P/(h*w)+&
            alpha*(-SIX*x*P*y**TWO/(h**THREE*w)+ &
            SIX*P*L*y**TWO/(h**THREE*w)+ &
            TWO*x**THREE*P/(h**THREE*w)- &
            SIX*P*L*x**TWO/(h**THREE*w)+ &
            HALF*(-TWO*P+alpha*P)*x/(w*alpha*h)))/shear_modulus

        gd3 = ZERO

        ivar = nint(self%get_option_value('ivar'))
        select case (ivar)
            case (1)
                val = gd1 

            case (2)
                val = gd2

            case (3)
                val = gd3

        end select

    end function compute
    !**********************************************************************************


end module fcn_cantilevered_beam
