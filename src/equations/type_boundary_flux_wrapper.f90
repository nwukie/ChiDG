module type_boundary_flux_wrapper
    use atype_boundary_flux, only: boundary_flux_t


    type, public :: boundary_flux_wrapper_t
        class(boundary_flux_t), allocatable  :: flux
    end type boundary_flux_wrapper_t


end module type_boundary_flux_wrapper
