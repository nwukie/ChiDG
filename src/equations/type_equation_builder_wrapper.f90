module type_equation_builder_wrapper
    use type_equation_builder,  only: equation_builder_t
    implicit none



    type, public :: equation_builder_wrapper_t

        class(equation_builder_t),  allocatable :: bld

    end type equation_builder_wrapper_t






end module type_equation_builder_wrapper
