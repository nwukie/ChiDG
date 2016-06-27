module fcn_sin
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, ONE, TWO, THREE, FIVE, PI
    use type_point,     only: point_t
    use type_function,  only: function_t
    implicit none





    !>  Sine function in y-coordinate.
    !!
    !!  \f$ f(t,\vec{x}) = mean + amp sin(freq*y)   \f$
    !!
    !!  Function parameters:
    !!
    !!  \f$ mean - \text{Mean value offset}   \\
    !!      amp  - \text{Amplitude of the sine function}        \\
    !!      freq - \text{Frequency of the sine function}   \f$
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-------------------------------------------------------------------------------------
    type, extends(function_t), public :: sin_f


    contains

        procedure   :: init
        procedure   :: compute

    end type sin_f
    !*************************************************************************************



contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(sin_f),  intent(inout)   :: self

        !
        ! Set function name
        !
        call self%add_name("A_p_Bsin(y)")


        !
        ! Set function options to default settings
        !
        call self%add_option('mean', 1._rk)
        call self%add_option('amplitude', 1._rk)
        call self%add_option('frequency', 2._rk * PI)



    end subroutine init
    !*************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !----------------------------------------------------------------------------------
    impure elemental function compute(self,time,coord) result(val)
        class(sin_f),  intent(inout)  :: self
        real(rk),           intent(in)  :: time
        type(point_t),      intent(in)  :: coord

        real(rk)                        :: val

        real(rk)    :: x,   y,   z, &
                       mean, amp, frequency

        !
        ! Get inputs and function parameters
        !
        x = coord%c1_
        y = coord%c2_
        z = coord%c3_

        mean      = self%get_option_value('mean')
        amp       = self%get_option_value('amplitude')
        frequency = self%get_option_value('frequency')


        !
        ! Compute function
        !
        val = mean + amp * sin(frequency * y)

    end function compute
    !***********************************************************************************


end module fcn_sin
