module type_flux_wrapper
    use atype_flux, only: flux_t


    type, public :: flux_wrapper_t
        class(flux_t), allocatable  :: item
    end type flux_wrapper_t


end module type_flux_wrapper
