module type_function_info
    use mod_kinds,  only: ik
    use type_seed,  only: seed_t


    !> Container for flux information.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !---------------------------------------------------------------------
    type, public :: function_info_t

        integer(ik)     :: type     ! Function type. ex: 'boundary_advective_flux', 'boundary_diffusive_flux', etc.
        integer(ik)     :: ifcn     ! Index of the function of the give type being computed.
        integer(ik)     :: idepend  ! Dependency index of a given element to the function. Chimera could be > 1.
        integer(ik)     :: idiff    ! Index of the direction being linearized. XI_MIN, ETA_MAX, DIAG, etc.
        type(seed_t)    :: seed     ! Indices of the element being linearized

    end type function_info_t
    !*********************************************************************









end module type_function_info
