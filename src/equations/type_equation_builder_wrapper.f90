module type_equation_builder_wrapper
    use type_equation_builder,  only: equation_builder_t
    implicit none



    !>  A wrapper for equation_builder_t's so they can be stored in an array.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/2/2016
    !!
    !!
    !---------------------------------------------------------------------------------
    type, public :: equation_builder_wrapper_t

        class(equation_builder_t),  allocatable :: bld

    end type equation_builder_wrapper_t
    !**********************************************************************************






end module type_equation_builder_wrapper
