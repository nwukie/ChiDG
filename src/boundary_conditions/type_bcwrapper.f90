module type_bcwrapper
    use atype_bc,   only: bc_t
    implicit none
    private


    type, public :: bcwrapper_t
        class(bc_t), allocatable    :: bc
    end type bcwrapper_t


end module type_bcwrapper
