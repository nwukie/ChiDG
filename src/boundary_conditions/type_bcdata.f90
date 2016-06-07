module type_bcdata
    use mod_kinds,          only: ik
    use type_bcwrapper,     only: bcwrapper_t
    use type_connectivity,  only: connectivity_t
    implicit none





    !>  This container is used to load information from a grid file and pass that back
    !!  to ChiDG so it can be initialized.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, public :: bcdata_t

        character(len=:),       allocatable :: domain_      !< Domain name the bcdata is associated with
        type(bcwrapper_t),      allocatable :: bcs(:)       !< Array of boundary conditions for each face
        !integer(ik),            allocatable :: bcface(:)    !< bc face index. XI_MIN, XI_MAX, ETA_MIN, etc.
        !integer(ik),            allocatable :: bc_connectivity(:,:)    !< bc face index. XI_MIN, XI_MAX, ETA_MIN, etc.
        type(connectivity_t),   allocatable :: bc_connectivity(:)    !< bc face index. XI_MIN, XI_MAX, ETA_MIN, etc.

    contains


    end type bcdata_t
    !******************************************************************************


















end module type_bcdata
