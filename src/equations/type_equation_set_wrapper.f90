module type_equation_set_wrapper
    use type_equation_set,  only: equation_set_t
    implicit none



    !>  A wrapper for equation_builder_t's so they can be stored in an array.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/2/2016
    !!
    !!
    !---------------------------------------------------------------------------------
    type, public :: equation_set_wrapper_t

        class(equation_set_t),  allocatable :: eqn

    end type equation_set_wrapper_t
    !**********************************************************************************






end module type_equation_set_wrapper
