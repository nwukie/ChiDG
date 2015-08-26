module type_bcwrapper
    use atype_bc,   only: bc_t
    implicit none
    private

    !>  Wrapper for storing a polymorphic boundary condition type bc_t
    !!      - This allows one to store an array of bcwrapper_t. A work around for storing an array
    !!        of polymorphic entities
    !!
    !!  @author Nathan A. Wukie
    !---------------------------------------------------
    type, public :: bcwrapper_t
        class(bc_t), allocatable    :: bc
    end type bcwrapper_t


end module type_bcwrapper
