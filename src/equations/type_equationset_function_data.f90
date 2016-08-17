module type_equationset_function_data
    implicit none



    !> Store information about the functions in an equation set
    !!
    !!  @author Nathan A. Wukie
    !!
    !-----------------------------------------------------------------
    type, public :: equationset_function_data_t

        integer :: nboundary_advective_flux = 0
        integer :: nboundary_diffusive_flux = 0

        integer :: nvolume_advective_flux   = 0
        integer :: nvolume_diffusive_flux   = 0

        integer :: nvolume_advective_source = 0
        integer :: nvolume_diffusive_source = 0

    end type equationset_function_data_t
    !******************************************************************




end module type_equationset_function_data
