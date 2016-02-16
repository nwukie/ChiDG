module type_function_wrapper
    use type_function,   only: function_t
    implicit none
    private

    !>  Wrapper for storing a polymorphic function type function_t
    !!      - This allows one to store an array of function_t. A work around for storing an array
    !!        of polymorphic entities
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/9/2016
    !!
    !-------------------------------------------------------------
    type, public :: function_wrapper_t

        class(function_t), allocatable    :: fcn

    end type function_wrapper_t
    !*************************************************************


end module type_function_wrapper
