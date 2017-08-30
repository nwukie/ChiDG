module fcn_turbine_bc_sine
#include<messenger.h>
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: ZERO,ONE,TWO,THREE,FOUR,FIVE,PI
    use type_point,     only: point_t
    use type_function,  only: function_t
    implicit none



    !>  Sine function in two coordinates
    !!
    !!  \f$ f(t, \vec{x}) = mean + amp*sin(freq_1*coord_1 + freq_2*coord_2 + phase) \f$
    !!
    !!  Function parameters:
    !!
    !!  mean    - Mean values offset
    !!  amp     - Amplitude of the sine function
    !!  freq_1  - Freqeuncy of the first coordinate of the sine function
    !!  freq_2  - Frequency of the second coordinate of the sine function
    !!  coord_1 - First coordinate of the sine function in [x,y,z,t]
    !!  coord_2 - Second coordinate of the sine function in [x,y,z,t]
    !!  phase   - Phase of the sine function
    !!
    !!  @author Mayank Sharma
    !!  @date   8/27/2017
    !!
    !-----------------------------------------------------------------------------------------------------------------
    type,   extends(function_t),    public      :: turbine_bc_sine_f
    

    contains

        procedure   :: init
        procedure   :: compute


    end type turbine_bc_sine_f
    !*****************************************************************************************************************



contains



    !>
    !!
    !!  @author Mayank Sharma
    !!  @date   8/27/2017
    !!
    !-----------------------------------------------------------------------------------------------------------------
    subroutine init(self)
        class(turbine_bc_sine_f),   intent(inout)       :: self

        !
        ! Set function name
        !
        call self%set_name("turbine_bc_sinusoid")

        !
        ! Set function to defualt settings
        !
        call self%add_option("mean", ONE)
        call self%add_option("amplitude", ONE)
        call self%add_option("frequency_1", TWO*PI)
        call self%add_option("frequency_2", TWO*PI)
        call self%add_option("coordinate_1", ONE)
        call self%add_option("coordinate_2", TWO)
        call self%add_option("phase", ZERO)


    end subroutine init
    !*****************************************************************************************************************



    !>
    !!
    !!  @author Mayank Sharma
    !!  @date   8/27/2017
    !!
    !-----------------------------------------------------------------------------------------------------------------
    impure elemental function compute(self,time,coord) result(val)
        class(turbine_bc_sine_f),   intent(inout)       :: self
        real(rk),                   intent(in)          :: time
        type(point_t),              intent(in)          :: coord

        real(rk)                                        :: val

        real(rk)        :: x,y,z,mean,amp,freq_1,freq_2,coord_1,coord_2, &
                           phase,compute_coord_1,compute_coord_2

        !
        ! Get inputs and function parameters
        !
        x = coord%c1_
        y = coord%c2_
        z = coord%c3_

        mean    = self%get_option_value("mean")
        amp     = self%get_option_value("amplitude")
        freq_1  = self%get_option_value("frequency_1")
        freq_2  = self%get_option_value("frequency_2")
        coord_1 = self%get_option_value("coordinate_1")
        coord_2 = self%get_option_value("coordinate_2")
        phase   = self%get_option_value("phase")

        !
        ! Compute coordinate 1
        !
        if (coord_1 == ONE) then
            compute_coord_1 = x
        else if (coord_1 == TWO) then
            compute_coord_1 = y
        else if (coord_1 == THREE) then
            compute_coord_1 = z
        else if (coord_1 == FOUR) then
            compute_coord_1 = time
        else
            compute_coord_1 = ZERO
        end if

        !
        ! Compute coordinate 2
        !
        if (coord_2 == ONE) then
            compute_coord_2 = x
        else if (coord_2 == TWO) then
            compute_coord_2 = y
        else if (coord_2 == THREE) then
            compute_coord_2 = z
        else if (coord_2 == FOUR) then
            compute_coord_2 = time
        else
            compute_coord_2 = ZERO
        end if
        
    !!  \f$ f(t, \vec{x}) = mean + amp*sin(freq_1*coord_1 + freq_2*coord_2 + phase) \f$
        val = mean + amp*sin(freq_1*compute_coord_1 + freq_2*compute_coord_2 + phase)

         
    end function compute
    !*****************************************************************************************************************




















end module fcn_turbine_bc_sine
