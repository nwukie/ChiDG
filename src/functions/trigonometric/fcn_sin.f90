module fcn_sine
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO,ONE,TWO,THREE,FIVE,PI
    use type_point,     only: point_t
    use type_function,  only: function_t
    implicit none




    !>  Sine function in a particular coordinate
    !!
    !!  \f$ f(t, \vec{x}) = mean + amp*sin(freq*coordinate + phase) \f$
    !!
    !!  Function parameters:
    !!
    !!  \f$ mean       - \text{Mean values offset}    \\
    !!      amp        - \text{Amplitude of the sine function}    \\
    !!      freq       - \text{Frequency of the sine function}    \\
    !!      coordinate - \text{Coordinate for the sine function \in [x,y,z,t]}  \\
    !!      phase      - \text{Phase shift of the sine function}    \f$
    !!
    !!  @author Mayank Sharma
    !!  @date   27/1/2016
    !!
    !------------------------------------------------------------------------------------------------------
    type, extends(function_t), public   :: sine_f


    contains

        procedure   :: init
        procedure   :: compute

    end type sine_f
    !******************************************************************************************************



contains



    !>
    !!
    !!  @author Mayank Sharma
    !!  @date   27/1/2017
    !!
    !------------------------------------------------------------------------------------------------------
    subroutine init(self)
        class(sine_f),  intent(inout)   :: self

        ! 
        ! Set function name
        !
        call self%set_name("sinusoid")

        !
        ! Set function to default settings
        !
        call self%add_option('mean', 1._rk)
        call self%add_option('amplitude', 1._rk)
        call self%add_option('frequency', 2._rk*PI)
        call self%add_option('coordinate',1._rk)
        call self%add_option('phase',0._rk)


    end subroutine init
    !******************************************************************************************************



    !>
    !!
    !!  @author Mayank Sharma
    !!  @date   27/1/2017
    !!
    !------------------------------------------------------------------------------------------------------
    impure elemental function compute(self,time,coord) result(val)
        class(sine_f),  intent(inout)      :: self
        real(rk),       intent(in)         :: time
        type(point_t),  intent(in)         :: coord

        real(rk)                           :: val

        real(rk)    :: x, y, z, mean, amp, frequency, coordinate, phase

        !
        ! Get inputs and function parameters
        !
        x = coord%c1_
        y = coord%c2_
        z = coord%c3_

        mean       = self%get_option_value('mean')
        amp        = self%get_option_value('amplitude')
        frequency  = self%get_option_value('frequency')
        coordinate = self%get_option_value('coordinate')
        phase      = self%get_option_value('phase')

        !
        ! Compute function
        !
        if (coordinate == 1.0_rk) then
            val = mean + amp*sin(frequency*x + phase)
        else if (coordinate == 2.0_rk) then
            val = mean + amp*sin(frequency*y + phase)
        else if (coordinate == 3.0_rk) then
            val = mean + amp*sin(frequency*z + phase)
        else if (coordinate == 4.0_rk) then
            val = mean + amp*sin(frequency*time + phase)
        else
            call chidg_signal(FATAL, "coordinate: invalid option value")
        end if


    end function compute
    !******************************************************************************************************




















end module fcn_sine
