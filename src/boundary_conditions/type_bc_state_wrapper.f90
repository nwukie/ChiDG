module type_bc_state_wrapper
    use type_bc_state,  only: bc_state_t
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !-----------------------------------------------------------------------
    type, public :: bc_state_wrapper_t

        class(bc_state_t),  allocatable :: state

    end type bc_state_wrapper_t
    !************************************************************************






end module type_bc_state_wrapper
