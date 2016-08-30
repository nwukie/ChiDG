module type_bc_operator_wrapper
    use type_bc_operator,  only: bc_operator_t
    implicit none


    !>  Wrapper for storing an array of dynamic equationset_t allocations.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   8/29/2016
    !!
    !!
    !----------------------------------------------------------------------------
    type, public :: bc_operator_wrapper_t

        class(bc_operator_t), allocatable  :: op

    end type bc_operator_wrapper_t
    !*****************************************************************************


end module type_bc_operator_wrapper
