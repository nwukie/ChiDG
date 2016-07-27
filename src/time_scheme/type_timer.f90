module type_timer
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: ZERO
    use mod_chidg_mpi,  only: GLOBAL_MASTER
    use mpi_f08,        only: MPI_WTime
    implicit none



    !>  Timer type for timing procedures
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !-----------------------------------------------------------
    type, public    :: timer_t

        real(rk)    :: start_time           !< Timer starting time
        real(rk)    :: stop_time            !< Timer ending time
        real(rk)    :: elapsed_time = ZERO  !< Timer elapsed time

        logical     :: started = .false.    !< Logical state indicating if the timer was started
        logical     :: stopped = .false.    !< Logical state indicating if the timer was stopped


    contains
        procedure   :: start                !< Start timer
        procedure   :: stop                 !< Stop timer
        procedure   :: elapsed              !< Compute elapsed time
        procedure   :: report               !< Print elapsed time
        procedure   :: reset                !< Reset timer

    end type timer_t
    !-----------------------------------------------------------



contains


    !> Start the timer
    !!
    !!  @author Nathan A. Wukie
    !!
    !---------------------------------------------------------------------------
    subroutine start(self)
        class(timer_t), intent(inout)  :: self

        if (.not. self%started) then
            !call cpu_time(self%start_time)
            self%start_time = MPI_WTime()
        else
            call chidg_signal(WARN,'type_timer: Timer was already started')
        end if

        self%started = .true.
        
    end subroutine start
    !***************************************************************************






    !> Stop the timer
    !!
    !!  @author Nathan A. Wukie
    !!
    !---------------------------------------------------------------------------
    subroutine stop(self)
        class(timer_t), intent(inout)  :: self

        if (self%started) then
            !call cpu_time(self%stop_time)
            self%stop_time = MPI_Wtime()
        else
            call chidg_signal(WARN,'type_timer: Timer was never started')
        end if


        self%elapsed_time = self%elapsed_time + real(self%stop_time - self%start_time,rk)

        self%stopped = .true.
        self%started = .false.

    end subroutine stop
    !***************************************************************************






    !> Return the elapsed time
    !!
    !!  @author Nathan A. Wukie
    !!
    !----------------------------------------------------------------------------
    function elapsed(self) result(elapsed_time)
        class(timer_t), intent(inout)  :: self

        real(rk)    :: elapsed_time
        
        if (self%stopped) then
            !elapsed_time = real(self%stop_time - self%start_time,rk)
            elapsed_time = self%elapsed_time
        else
            call chidg_signal(WARN,'type_timer: Timer was not stopped')
            elapsed_time = 123456789._rk
        end if

    end function elapsed
    !****************************************************************************






    !> Procedure for reporting timer information
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  string  Optional character string for reporting contextual information about the timing
    !!
    !------------------------------------------------------------------------------------------------------
    subroutine report(self,messg)
        class(timer_t), intent(inout)           :: self
        character(*),   intent(in),  optional   :: messg

        real(rk)    :: elapsed


        ! Compute the elapsed time
        elapsed =  self%elapsed()


        ! If a message was included, print and append timing. Else, just print the timing
        if (present(messg)) then
            call write_line(messg, ' ', elapsed, delimiter='', io_proc=GLOBAL_MASTER)
        else
            call write_line(elapsed, io_proc=GLOBAL_MASTER)
        end if


    end subroutine report
    !*******************************************************************************************************




    !> Procedure for reseting the clock
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !-------------------------------------------------------------------
    subroutine reset(self)
        class(timer_t), intent(inout)       :: self

        self%start_time   = 0.
        self%stop_time    = 0.
        self%elapsed_time = 0.

        self%started      = .false.
        self%stopped      = .false.

    end subroutine reset
    !********************************************************************************************************




end module type_timer
