module mod_timescheme
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use type_timescheme,    only: timescheme_t
    use type_dict,          only: dict_t




    ! Import solverdata types
!    use forward_euler,                  only: forward_euler_t
!    use backward_euler,                 only: backward_euler_t
!    use backward_euler_subiteration,    only: backward_euler_subiteration_t
    use quasi_newton,                   only: quasi_newton_t
    use newton,                         only: newton_t
    implicit none



    ! Instantiate solver types for sourcing
!    type(forward_euler_t)               :: FORWARD_EULER_INSTANCE
!    type(backward_euler_t)              :: BACKWARD_EULER_INSTANCE
!    type(backward_euler_subiteration_t) :: BACKWARD_EULER_SUBITERATION_INSTANCE
    type(quasi_newton_t)                :: QUASI_NEWTON_INSTANCE
    type(newton_t)                      :: NEWTON_INSTANCE



    logical :: initialized = .false.



contains





    subroutine create_timescheme(time_string,instance,options)
        character(*),                       intent(in)      :: time_string
        class(timescheme_t), allocatable,   intent(inout)   :: instance
        type(dict_t), optional,             intent(inout)   :: options

        select case (trim(time_string))
!            case ('ForwardEuler','forwardeuler','Forward_Euler','forward_euler','FE','fe')
!                allocate(instance, source=FORWARD_EULER_INSTANCE)
!
!
!            case ('BackwardEuler','backwardeuler','Backward_Euler','backward_euler','BE','be')
!                allocate(instance, source=BACKWARD_EULER_INSTANCE)
!
!            case ('BackwardEulerSub','backwardeulersub','Backward_Euler_Sub','backward_euler_sub','backward_euler_subiteration','BES','bes')
!                allocate(instance, source=BACKWARD_EULER_SUBITERATION_INSTANCE)
!
            case ('quasi_newton','Quasi_Newton','quasinewton','QuasiNewton')
                allocate(instance, source=QUASI_NEWTON_INSTANCE)

            case ('newton','Newton')
                allocate(instance, source=NEWTON_INSTANCE)


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


    end subroutine








end module mod_timescheme
