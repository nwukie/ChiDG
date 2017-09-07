module type_tutorial_step_wrapper
    use type_tutorial_step, only: tutorial_step_t
    implicit none

    type, public :: tutorial_step_wrapper_t
        class(tutorial_step_t), allocatable :: item
    end type tutorial_step_wrapper_t

end module type_tutorial_step_wrapper
