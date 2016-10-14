module mod_time_scheme
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use type_time_scheme,   only: time_scheme_t
    use type_dict,          only: dict_t




    ! Import solverdata types
    use steady,                     only: steady_t
    use forward_euler,              only: forward_euler_t
    implicit none



    ! Instantiate solver types for sourcing
    type(steady_t)                      :: STEADY
    type(forward_euler_t)               :: FORWARD_EULER



    logical :: initialized = .false.



contains




    !> Create a concrete timescheme
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/15/2016
    !!
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------------------------------
    subroutine create_time_scheme(time_string,instance,options)
        character(*),                       intent(in)      :: time_string
        class(time_scheme_t), allocatable,  intent(inout)   :: instance
        type(dict_t), optional,             intent(inout)   :: options

        character(len=:), allocatable   :: user_msg, dev_msg



        select case (trim(time_string))

            case ('steady','Steady','STEADY')
                allocate(instance, source=STEADY)

            case ('forward_euler','Forward_Euler','FORWARD_EULER','forward euler')
                allocate(instance, source=FORWARD_EULER)

            case default
                user_msg = "We can't seem to find a time integrator that matches the input string. &
                            Maybe check that the time integrator string in the input file or driver &
                            script is valid."
                dev_msg = "Check that the time integrator is registered properly in create_time_scheme."
                call chidg_signal_two(OOPS, user_msg, trim(time_string), dev_msg=dev_msg)
        end select





        !
        ! Call options initialization if present
        !
        if (present(options)) then
            call instance%set(options)
        end if



        !
        ! Make sure the solver was allocated
        !
        if (.not. allocated(instance)) call chidg_signal(FATAL,"create_time_scheme: solver was not allocated. Check that the desired solver was registered and instantiated in the mod_time_scheme module")


    end subroutine create_time_scheme
    !***********************************************************************************************************************








end module mod_time_scheme
