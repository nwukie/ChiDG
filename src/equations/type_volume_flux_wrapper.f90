module type_volume_flux_wrapper
    use atype_volume_flux, only: volume_flux_t


    type, public :: volume_flux_wrapper_t
        class(volume_flux_t), allocatable  :: flux
    end type volume_flux_wrapper_t


end module type_volume_flux_wrapper
