module type_equationset_wrapper
    use type_equationset,   only: equationset_t
    implicit none


    !>  Wrapper for storing an array of dynamic equationset_t allocations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !----------------------------------------------------------------------------
    type, public :: equationset_wrapper_t

        class(equationset_t), allocatable   :: item

    end type equationset_wrapper_t
    !*****************************************************************************


end module type_equationset_wrapper
