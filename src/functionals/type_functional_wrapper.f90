module type_functional_wrapper
    use type_evaluator,    only:evaluator_t
    implicit none


    !>
    !!
    !!  @author Matteo Ugolotti
    !!  @data   05/10/2017
    !!
    !---------------------------------------------------------------------
    type,public :: functional_wrapper_t
        class(evaluator_t),    allocatable :: func
    end type functional_wrapper_t
    !*********************************************************************


end module type_functional_wrapper
