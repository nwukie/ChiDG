module mod_time_scheme
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use type_time_scheme,   only: time_scheme_t
    use type_dict,          only: dict_t




    ! Import solverdata types
    use steady,                         only: steady_t
    implicit none



    ! Instantiate solver types for sourcing
    type(steady_t)                      :: STEADY



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





        select case (trim(time_string))

            case ('steady','Steady')
                allocate(instance, source=STEADY)


            case default
                call chidg_signal(FATAL,'create_time_scheme -- solver string not recognized')
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
