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

        ! Instantiate time_integrators
        type(steady_t)                      :: steady
        type(forward_euler_t)               :: forward_euler
        type(backward_euler_t)              :: backward_euler
        type(harmonic_balance_t)            :: harmonic_balance
        type(explicit_runge_kutta_t)        :: explicit_runge_kutta
        type(DIRK_t)                        :: DIRK
        type(DIRK_coupled_oscillator_t)     :: DIRK_coupled_oscillator


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

        ! Find equation set in 'available_equations' vector
        bindex = self%time_integrators%index_by_name(string)

        ! Check time integrator was found in factory
        user_msg = "We can't seem to find a time integrator that matches the string that &
                    was passed into the time integrator factory. Maybe check that the &
                    time integrator strings that were set for the domains are all valid."
        dev_msg = "Check that the time integrator is registered properly in the &
                   time integrator factory: &
                   src/time_integrator/mod_time_integrators.f90. When a time integrator &
                   is defined, is needs to be registered in the factory by &
                   calling 'call time_integrator_factory%register(time_integrator)', &
                   where 'time_integrator' is the time_integrator instance to be registered. &
                   This could be done in 'mod_time_integrators.register_time_integrators'. &
                   The 'register_time_integrators' routine gets called on startup and &
                   loads all the default time integrators into the factory so the library &
                   knows what time integrators it can build."
        if (bindex == 0) call chidg_signal_two(OOPS,user_msg,trim(string),dev_msg=dev_msg)

        ! Get equation set builder
        allocate(time_integrator, source=self%time_integrators%at(bindex), stat=ierr)
        if (ierr /= 0) call AllocationError

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



end module mod_time_integrators
