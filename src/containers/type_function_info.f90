module type_function_info
    use mod_kinds,  only: ik


    !> Container for flux information.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !---------------------------------------------------------------------
    type, public :: function_info_t

        integer(ik)             :: type     !< Function type. ex: 'boundary_advective_flux', 'boundary_diffusive_flux', etc.
        integer(ik)             :: ifcn     !< Index of the function of the give type being computed.
        integer(ik)             :: idonor   !< Index of the donor to the face. Chimera could be > 1.
        integer(ik)             :: iblk     !< Index of the block being linearized.

    end type function_info_t
    !*********************************************************************









end module type_function_info
