module mod_time
#include<messenger.h>
    use type_time_manager,  only: time_manager_t

    implicit none

    type(time_manager_t)    :: time_manager



contains



    !>
    !!
    !!
    !!
    !!
    !--------------------------------------------------------
    subroutine initialize_time_manager()

        call time_manager%init()


    end subroutine initialize_time_manager
    !********************************************************



















end module mod_time
