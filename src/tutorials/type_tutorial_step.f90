module type_tutorial_step
#include <messenger.h>
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/17/2017
    !!
    !!
    !-----------------------------------------------------------------------------------
    type, public :: tutorial_step_t

        character(:),               allocatable :: title
        class(tutorial_action_t),   allocatable :: action

    contains

        procedure   :: set_title
        procedure   :: get_title
        procedure(deferred)   :: execute

    end type tutorial_step_t
    !***********************************************************************************




contains


    !>  Set the title of the step.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/17/2017
    !!
    !-----------------------------------------------------------------------------------
    subroutine set_title(self,title)
        class(tutorial_step_t), intent(inout)   :: self
        character(*),           intent(in)      :: title

        self%title = title

    end subroutine set_title
    !***********************************************************************************




    !>  Return the title of the step.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/17/2017
    !!
    !-----------------------------------------------------------------------------------
    function get_title(self) result(title)
        class(tutorial_step_t), intent(in)  :: self

        character(:),   allocatable :: title

        title = self%title

    end function get_title
    !***********************************************************************************



    !>  Execute the step procedures.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/17/2017
    !!
    !-----------------------------------------------------------------------------------
    subroutine execute(self)
        class(tutorial_step_t), intent(inout)   :: self



    end subroutine execute
    !***********************************************************************************



end module type_tutorial_step
