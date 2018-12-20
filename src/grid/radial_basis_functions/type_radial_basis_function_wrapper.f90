module type_radial_basis_function_wrapper
    use type_radial_basis_function,   only: radial_basis_function_t
    implicit none
    private

    !>  Wrapper for storing a polymorphic function type radial_basis_function_t
    !!      - This allows one to store an array of radial_basis_function_t. A work around for storing an array
    !!        of polymorphic entities
    !!
    !!  @author Eric Wolf
    !!  @date   3/15/2017
    !!
    !-------------------------------------------------------------
    type, public :: radial_basis_function_wrapper_t

        class(radial_basis_function_t), allocatable    :: rbf

    end type radial_basis_function_wrapper_t
    !*************************************************************


end module type_radial_basis_function_wrapper
