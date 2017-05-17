module type_tutorial
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use type_tutorial_step, only: tutorial_step_t
    implicit none



    !>  An object for driving tutorials.
    !!
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------
    type, abstract, public :: tutorial_t

        character(:),   allocatable         :: title
        character(:),   allocatable         :: tags
        character(:),   allocatable         :: files

        type(tutorial_step_t),  allocatable ::  step(:)

    contains

        procedure(self_interface),  deferred    :: print_configuration
        procedure(self_interface),  deferred    :: initialize

        procedure                               :: set_title
        procedure                               :: get_title
        procedure                               :: set_tags
        procedure                               :: get_tags
        procedure                               :: set_files
        procedure                               :: get_files
        procedure                               :: add_step
        procedure                               :: new_step
        procedure                               :: nsteps
        procedure                               :: run
        procedure                               :: stepper

    end type tutorial_t
    !***********************************************************************




    abstract interface
        subroutine self_interface(self)
            import tutorial_t
            class(tutorial_t),  intent(inout)   :: self
        end subroutine
    end interface




contains


    !>
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_title(self,title)
        class(tutorial_t),  intent(inout)   :: self
        character(*),       intent(in)      :: title

        self%title = trim(title)

    end subroutine set_title
    !****************************************************************************************



    !>
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_title(self) result(title)
        class(tutorial_t),  intent(in)  :: self

        character(:),   allocatable :: title

        title = self%title

    end function get_title
    !****************************************************************************************



    !>
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_tags(self,tags)
        class(tutorial_t),  intent(inout)   :: self
        character(*),       intent(in)      :: tags

        self%tags = trim(tags)

    end subroutine set_tags
    !****************************************************************************************


    !>
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_tags(self) result(tags)
        class(tutorial_t),  intent(in)  :: self

        character(:),   allocatable :: tags

        tags = self%tags

    end function get_tags
    !****************************************************************************************



    !>
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine set_files(self,files)
        class(tutorial_t),  intent(inout)   :: self
        character(*),       intent(in)      :: files

        self%files = trim(files)

    end subroutine set_files
    !****************************************************************************************



    !>
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_files(self) result(files)
        class(tutorial_t),  intent(in)  :: self

        character(:),   allocatable :: files

        files = self%files

    end function get_files
    !****************************************************************************************



    !>  Add a step to the tutorial.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/17/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine add_step(self,step)
        class(tutorial_t),      intent(inout)   :: self
        class(tutorial_step_t), intent(in)      :: step

        integer(ik) :: step_ID

        step_ID = self%new_step()
        self%step(step_ID) = step 


    end subroutine add_step
    !****************************************************************************************




    !>  Allocate storage for new step and return ID.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/17/2017
    !!
    !----------------------------------------------------------------------------------------
    function new_step(self) result(step_ID)
        class(tutorial_t),  intent(inout)   :: self

        type(tutorial_step_t), allocatable  :: temp_steps(:)
        integer(ik)                         :: step_ID, ierr

        !
        ! Allocate number of steps
        !
        allocate(temp_steps(self%nsteps() + 1), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Copy any previously allocated steps to new array
        !
        if (self%nsteps() > 0) then
            temp_steps(1:self%nsteps()) = self%step(1:self%nsteps())
        end if

        
        !
        ! Set ID of new step to last place.
        !
        step_ID = size(temp_steps)


        !
        ! Attach extended allocation to self%step
        !
        call move_alloc(temp_steps,self%step)


    end function new_step
    !****************************************************************************************



    !>  Return the number of steps in the tutorial.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/17/2016
    !!
    !----------------------------------------------------------------------------------------
    function nsteps(self) result(nsteps_)
        class(tutorial_t),  intent(in)  :: self

        integer(ik) :: nsteps_

        if (allocated(self%step)) then
            nsteps_ = size(self%step)
        else
            nsteps_ = 0
        end if

    end function nsteps
    !****************************************************************************************




    !>  Run the tutorial.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/17/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine run(self)
        class(tutorial_t),  intent(inout)   :: self

        integer(ik)     :: step
        character(1024) :: user_input
        logical         :: done, read_input


        step = 1
        done = .false.
        do

            ! Header
            if (.not. done) then
                call execute_command_line('clear')
                call write_line('Title:',self%get_title(), bold=.true.)
                call write_line('Tags:',self%get_tags(), bold=.true.)
                call write_line('Files:',self%get_files(), bold=.true.)
                call self%print_configuration()
                call write_line('Step ', step,': ',self%step(step)%get_title(), color='blue', bold=.true.)
            end if



            ! Stepper - execute step
            call self%stepper(step,done)






            ! Detect last step and exit
            if (done) then
                call execute_command_line('clear')
                exit
            end if



        
            ! Wait for user to indicate next step
            read_input = .true.
            do while (read_input)

                read(*,'(A1024)') user_input

                ! Go forward/backward
                if (trim(user_input) == 'f') then
                    step = step + 1
                    read_input = .false.
                else if (trim(user_input) == 'b') then
                    step = step - 1
                    read_input = .false.
                else if (trim(user_input) == 'q') then
                    step = 0
                    read_input = .false.
                end if

            end do


            !
            ! Check if done
            !
            done =  (step==0) .or. (step > self%nsteps())


        end do

    end subroutine run
    !****************************************************************************************







    !>
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine stepper(self,step,done)
        class(tutorial_t),  intent(inout)   :: self
        integer(ik),        intent(in)      :: step
        logical,            intent(inout)   :: done


        if (.not. done) then
            call self%step(step)%execute()
        end if


    end subroutine stepper
    !****************************************************************************************




end module type_tutorial
