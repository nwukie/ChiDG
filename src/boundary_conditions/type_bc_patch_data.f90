module type_bc_patch_data
    use mod_kinds,                  only: ik
    use type_bcvector,              only: bcvector_t
    use type_svector,               only: svector_t
    use type_boundary_connectivity, only: boundary_connectivity_t
    implicit none





    !>  This container is used to load information from a grid file and pass that back
    !!  to ChiDG so it can be initialized.
    !!
    !!  Contains:
    !!      name of block the data is associated with
    !!      bc_connectivity: one for each boundary of the block. (6 total)
    !!      bc_group_name: vector containing name of bc_group each boundary 
    !!                     of the block is associated with.
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
        character(:),                   allocatable :: domain_name          ! Domain name the bcdata is associated with.
        type(boundary_connectivity_t),  allocatable :: bc_connectivity(:)   ! Face connectivities for faces defining the patch.
        type(svector_t)                             :: bc_group_name        ! Boundary State Group each face patch is associated with.

    end type bc_patch_data_t
    !******************************************************************************


















end module type_bc_patch_data
