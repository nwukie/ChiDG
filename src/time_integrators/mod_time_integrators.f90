module mod_time_integrators
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use type_time_integrator,   only: time_integrator_t
    use type_dict,              only: dict_t
    use type_tivector,          only: tivector_t
    use mod_time,               only: time_manager_global


    ! Import solverdata types
    use type_steady,                    only: steady_t
    use type_forward_euler,             only: forward_euler_t
    use type_backward_euler,            only: backward_euler_t
    use type_harmonic_balance,          only: harmonic_balance_t 
    use type_explicit_runge_kutta,      only: explicit_runge_kutta_t
    use type_DIRK,                      only: DIRK_t
    use type_DIRK_coupled_oscillator,   only: DIRK_coupled_oscillator_t
    implicit none



!    ! Instantiate solver types for sourcing
!    type(steady_t)                      :: STEADY
!    type(forward_euler_t)               :: FORWARD_EULER
!    type(backward_euler_t)              :: BACKWARD_EULER
!    type(harmonic_balance_t)            :: HB 
!    type(explicit_runge_kutta_t)        :: EXPLICIT_RK
!    type(DIRK_t)                        :: DIRK
!    type(DIRK_coupled_oscillator_t)     :: DIRK_coupled_oscillator





    !>  Factory object for managing creation of time_integrator_t.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/11/2017
    !!
    !--------------------------------------------------------------------------------------
    type, public :: time_integrator_factory_t

        type(tivector_t) :: time_integrators

    contains
        
        procedure   :: register => register_time_integrator
        procedure   :: produce  => produce_time_integrator
        procedure   :: list     => list_time_integrators
        procedure   :: has      => has_time_integrator

    end type time_integrator_factory_t
    !**************************************************************************************

    type(time_integrator_factory_t),   target   :: time_integrator_factory
    logical                                     :: initialized = .false.



contains


    !>  Register time_integrators in the time_integrator_factory. This is called 
    !!  from chidg%start_up('core')
    !!
    !!  This allows the available time_integrators to be queried in the same way 
    !!  that they are registered for allocation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/11/2017
    !!
    !-----------------------------------------------------------------------------------
    subroutine register_time_integrators()

        !
        ! Instantiate time_integrators
        !
        type(steady_t)                      :: steady
        type(forward_euler_t)               :: forward_euler
        type(backward_euler_t)              :: backward_euler
        type(harmonic_balance_t)            :: harmonic_balance
        type(explicit_runge_kutta_t)        :: explicit_runge_kutta
        type(DIRK_t)                        :: DIRK
        type(DIRK_coupled_oscillator_t)     :: DIRK_coupled_oscillator


        !
        ! Register if needed
        !
        if ( .not. initialized ) then

            ! Register in global vector
            call time_integrator_factory%register(steady)
            call time_integrator_factory%register(forward_euler)
            call time_integrator_factory%register(backward_euler)
            call time_integrator_factory%register(harmonic_balance)
            call time_integrator_factory%register(explicit_runge_kutta)
            call time_integrator_factory%register(DIRK)
            call time_integrator_factory%register(DIRK_coupled_oscillator)

            ! Confirm initialization
            initialized = .true.

        end if

    end subroutine register_time_integrators
    !*************************************************************************************







    !>  Register a time_integrator with the factory.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/11/2017
    !!
    !--------------------------------------------------------------------------------------
    subroutine register_time_integrator(self,time_integrator)
        class(time_integrator_factory_t),   intent(inout)   :: self
        class(time_integrator_t),           intent(inout)   :: time_integrator

        call time_integrator%init()

        call self%time_integrators%push_back(time_integrator)

    end subroutine register_time_integrator
    !***************************************************************************************



    !>  time_integrator_t factory:
    !!      - procedure for allocating a concrete instance of an time_integrator_t
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/11/2017
    !!
    !!  @param[in] string    Character string for the equation set name
    !!
    !-------------------------------------------------------------------------------------
    function produce_time_integrator(self,string) result(time_integrator)
        class(time_integrator_factory_t),   intent(inout)   :: self
        character(*),                       intent(in)      :: string

        character(:),               allocatable :: user_msg, dev_msg
        integer                                 :: ierr, bindex
        class(time_integrator_t),   allocatable :: time_integrator


        !
        ! Find equation set in 'available_equations' vector
        !
        bindex = self%time_integrators%index_by_name(string)


        !
        ! Check time integrator was found in factory
        !
        user_msg = "We can't seem to find an equation set that matches the string that &
                    was passed into the equation set factory. Maybe check that the &
                    equation set strings that were set for the domains are all valid."
        dev_msg = "Check that the equation set builder is registered properly in the &
                   equation set builder factory: &
                   src/time_integrator/mod_time_integrators.f90. When an equation &
                   set builder is defined, is needs to be registered in the factory by &
                   calling 'call time_integrator_factory%register(time_integrator)', &
                   where 'builder' is the object that knows how to build the equation set. &
                   This could be done in 'mod_equations.register_equation_builders'. &
                   The 'register_equation_builders' routine gets called on startup and &
                   loads all the default equation builders into the factory so the library &
                   knows what equation sets it can build."
        if (bindex == 0) call chidg_signal_two(OOPS,user_msg,trim(string),dev_msg=dev_msg)


        !
        ! Get equation set builder
        !
        allocate(time_integrator, source=self%time_integrators%at(bindex), stat=ierr)
        if (ierr /= 0) call AllocationError
        !time_integrator = self%time_integrators%at(bindex)


    end function produce_time_integrator
    !*************************************************************************************





    !>  This is really a utilitity for 'chidg edit' to dynamically list the avalable 
    !!  equation sets.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/11/2017
    !!
    !--------------------------------------------------------------------------------------
    subroutine list_time_integrators(self)
        class(time_integrator_factory_t),   intent(in)  :: self

        integer(ik)                 :: item
        character(:),   allocatable :: tname

        do item = 1,self%time_integrators%size()
            tname = self%time_integrators%data(item)%item%get_name()
            call write_line(trim(tname))
        end do ! item

    end subroutine list_time_integrators
    !**************************************************************************************





    !>  Check if a given equation builder is registered in the factory.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2017
    !!
    !--------------------------------------------------------------------------------------
    function has_time_integrator(self,string) result(integrator_status)
        class(time_integrator_factory_t),   intent(in)  :: self
        character(*),                       intent(in)  :: string

        integer(ik)                 :: ind
        character(:),   allocatable :: tname
        logical                     :: integrator_status
        
        integrator_status = .false.
        do ind = 1,self%time_integrators%size()
            tname = self%time_integrators%data(ind)%item%get_name()
            integrator_status = (trim(string) == trim(tname))
            if (integrator_status) exit
        end do ! item

    end function has_time_integrator
    !**************************************************************************************










!    !>  Create a concrete time integrator
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   3/15/2016
!    !!
!    !!
!    !!
!    !------------------------------------------------------------------------------------------
!    subroutine create_time_integrator(time_string,instance)
!        character(*),                               intent(in)      :: time_string
!        class(time_integrator_t),   allocatable,    intent(inout)   :: instance
!
!        character(:),   allocatable :: user_msg, dev_msg
!
!
!
!        select case (trim(time_string))
!
!            case ('steady','Steady','STEADY')
!                allocate(instance, source=STEADY)
!
!            case ('forward_euler','Forward_Euler','FORWARD_EULER','forward euler', 'Forward Euler')
!                allocate(instance, source=FORWARD_EULER)
!
!            case ('backward_euler', 'Backward_Euler', 'BACKWARD_EULER', 'backward euler', 'Backward Euler', 'BACKWARD EULER')
!                allocate(instance, source=BACKWARD_EULER)
!            
!            case ('Harmonic Balance', 'Harmonic_Balance', 'harmonic balance', 'harmonic_balance', 'HB')
!                allocate(instance, source=HB)
!            
!            case ('Second Order Runge_Kutta', 'Explict Midpoint', 'Second Order RK', 'Modified Euler', 'Second Order Ralston Method', 'Third Order Runge-Kutta', 'Third Order Kutta', 'Third Order RK', &
!                   'Runge-Kutta Method', 'Fourth Runge-Kutta Method', 'Fourth Order RK Method', 'RK4', 'Three-Eighth Rule', 'Fourth Order Kutta') ! this probably needs to be split up in several RK schemes
!                allocate(instance, source=EXPLICIT_RK)
!
!            case ('DIRK')
!                allocate(instance, source=DIRK)
!
!            case ('DIRK_coupled_oscillator')
!                allocate(instance, source=DIRK_coupled_oscillator)
!
!            case default
!                user_msg = "We can't seem to find a time integrator that matches the input &
!                            string. Maybe check that the time integrator string in the input &
!                            file or driver script is valid."
!                dev_msg = "Check that the time integrator is registered properly in &
!                           create_time_integrator."
!                call chidg_signal_two(OOPS, user_msg, trim(time_string), dev_msg=dev_msg)
!        end select
!
!
!
!
!        !
!        ! Make sure the solver was allocated
!        !
!        user_msg = "create_time_integrator: solver was not allocated. Check that the desired &
!                                            solver was registered and instantiated in the mod_time_integrator module"
!        if (.not. allocated(instance)) call chidg_signal(FATAL,user_msg)
!
!
!    end subroutine create_time_integrator
!    !****************************************************************************************








end module mod_time_integrators
