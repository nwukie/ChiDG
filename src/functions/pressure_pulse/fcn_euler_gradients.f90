module fcn_euler_gradients
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
    type, extends(function_t), public :: euler_gradients_f

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

    end type euler_gradients_f
    !********************************************************************************



contains




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(euler_gradients_f),  intent(inout)  :: self

        !
        ! Set function name
        !
        call self%set_name("euler_gradients")

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
        class(euler_gradients_f),     intent(inout)  :: self
        real(rk),                       intent(in)  :: time
        type(point_t),                  intent(in)  :: coord

        real(rk)                                    :: val
        integer(ik)                                 :: ivar

        real(rk)    :: x, y, z

        x   = coord%c1_
        y   = coord%c2_
        z   = coord%c3_

        ivar = nint(self%get_option_value('ivar'))
        select case (ivar)
            ! RHO
            case (1)
                val = x

            ! RHO-U
            case (2)
                val = 5.0*x*x

            ! RHO-V
            case (3)
                val = ZERO

            ! RHO-W
            case (4)
                val = ZERO

            ! RHO-E
            case (5)
                val = 10.0*x

        end select

    end function compute
    !**********************************************************************************


end module fcn_euler_gradients
