module mod_tutorials
#include <messenger.h>
!    use mod_tutorial_smoothbump,    only: tutorial_driver_smoothbump
    use mod_tutorial_smoothbump_new,    only: smooth_bump_tutorial_t
    implicit none




contains


    !>
    !!
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------
    subroutine tutorial_driver(tutorial)
        character(*),   intent(in)  :: tutorial

        type(smooth_bump_tutorial_t)    :: smooth_bump_tutorial


        select case (trim(tutorial))
            case ('smooth_bump')
                !call tutorial_driver_smoothbump()
                call smooth_bump_tutorial%initialize()
                call smooth_bump_tutorial%run()

            case default
                call chidg_signal(FATAL,"tutorial_driver: the available tutorials are 'smooth_bump'.")

        end select


    end subroutine tutorial_driver
    !************************************************************************




end module mod_tutorials
