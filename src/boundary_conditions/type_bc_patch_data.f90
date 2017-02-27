module type_bc_patch_data
    use mod_kinds,                  only: ik
    use type_bcvector,              only: bcvector_t
    use type_svector,               only: svector_t
    use type_boundary_connectivity, only: boundary_connectivity_t
    implicit none





    !>  This container is used to load information from a grid file and pass that back
    !!  to ChiDG so it can be initialized.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/31/2016
    !!  @note   Modified to use bcvector_t for holding vector of bc_operators for each face
    !!
    !------------------------------------------------------------------------------
    type, public :: bc_patch_data_t

        ! Boundary condition patch information
        character(:),                   allocatable :: domain_              ! Domain name the bcdata is associated with.
        type(boundary_connectivity_t),  allocatable :: bc_connectivity(:)   ! Face connectivities for faces defining the patch.
        type(svector_t)                             :: bc_group             ! Boundary State Group the patch is associated with.

    end type bc_patch_data_t
    !******************************************************************************


















end module type_bc_patch_data
