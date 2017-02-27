module mod_time_integrators
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use type_time_integrator,   only: time_integrator_t
    use type_dict,              only: dict_t
    use mod_time,               only: time_manager_global


    ! Import solverdata types
    use type_steady,               only: steady_t
    use type_forward_euler,        only: forward_euler_t
!    use type_backward_euler,    only: backward_euler_t
    use type_harmonic_balance,     only: harmonic_balance_t !not in yet
    use type_explicit_runge_kutta, only: explicit_runge_kutta_t
    implicit none



    ! Instantiate solver types for sourcing
    type(steady_t)                      :: STEADY
    type(forward_euler_t)               :: FORWARD_EULER
!    type(backward_euler_t)              :: BACKWARD_EULER
    type(harmonic_balance_t)            :: HB !not in yet
    type(explicit_runge_kutta_t)        :: EXPLICIT_RK

    logical :: initialized = .false.



contains




    !>  Create a concrete time integrator
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/15/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine create_time_integrator(time_string,instance)
        character(*),                               intent(in)      :: time_string
        class(time_integrator_t),   allocatable,    intent(inout)   :: instance

        character(:),   allocatable :: user_msg, dev_msg



        select case (trim(time_string))

            case ('steady','Steady','STEADY')
                allocate(instance, source=STEADY)

            case ('forward_euler','Forward_Euler','FORWARD_EULER','forward euler', 'Forward Euler')
                allocate(instance, source=FORWARD_EULER)
!
!            case ('backward_euler', 'Backward_Euler', 'BACKWARD_EULER', 'backward euler', 'Backward Euler', 'BACKWARD EULER')
!                allocate(instance, source=BACKWARD_EULER)
            
            case ('Harmonic Balance', 'Harmonic_Balance', 'harmonic balance', 'harmonic_balance', 'HB')
                allocate(instance, source=HB)
            
            case ('Second Order Runge_Kutta', 'Explict Midpoint', 'Second Order RK', 'Modified Euler', 'Second Order Ralston Method', 'Third Order Runge-Kutta', 'Third Order Kutta', 'Third Order RK', &
                   'Runge-Kutta Method', 'Fourth Runge-Kutta Method', 'Fourth Order RK Method', 'RK4', 'Three-Eighth Rule', 'Fourth Order Kutta') ! this probably needs to be split up in several RK schemes
                allocate(instance, source=EXPLICIT_RK)



            case default
                user_msg = "We can't seem to find a time integrator that matches the input &
                            string. Maybe check that the time integrator string in the input &
                            file or driver script is valid."
                dev_msg = "Check that the time integrator is registered properly in &
                           create_time_integrator."
                call chidg_signal_two(OOPS, user_msg, trim(time_string), dev_msg=dev_msg)
        end select





        !
        ! Call time_integrator initialization
        !
!        call instance%init() ! The time_integrator is now initialized in chidg%init() case 'finalized'
!        call instance%time_manager%init() ! the time_manager is now initialized in chidg%start_up case 'core'

!        if (present(options)) then
!            call instance%set(options)
!        end if



        !
        ! Make sure the solver was allocated
        !
        user_msg = "create_time_integrator: solver was not allocated. Check that the desired &
                                            solver was registered and instantiated in the mod_time_integrator module"
        if (.not. allocated(instance)) call chidg_signal(FATAL,user_msg)


    end subroutine create_time_integrator
    !****************************************************************************************








end module mod_time_integrators
