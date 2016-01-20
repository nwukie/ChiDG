module type_flux_indices
    use mod_kinds,  only: ik


    !> Container for flux information.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !---------------------------------------------------------------------
    type, public :: flux_indices_t

        integer(ik)     :: iflux    !< Index of the flux being computed.
        integer(ik)     :: idonor   !< Index of the donor to the face. Chimera could be > 1.
        integer(ik)     :: iblk     !< Index of the block being linearized.

    end type flux_indices_t
    !*********************************************************************









end module type_flux_indices
