module type_operator_wrapper
    use type_operator,  only: operator_t
    implicit none


    !>  Wrapper for storing an array of dynamic equationset_t allocations.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !----------------------------------------------------------------------------
    type, public :: operator_wrapper_t

        class(operator_t), allocatable  :: op

    end type operator_wrapper_t
    !*****************************************************************************


end module type_operator_wrapper
