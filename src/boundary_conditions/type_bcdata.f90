module type_bcdata
    use mod_kinds,                  only: ik
    use type_bcvector,              only: bcvector_t
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
    type, public :: bcdata_t

        character(len=:),               allocatable :: domain_              !< Domain name the bcdata is associated with
        type(bcvector_t),               allocatable :: bcs(:)               !< Vector of boundary condition operators for each face
        type(boundary_connectivity_t),  allocatable :: bc_connectivity(:)   !< 

    end type bcdata_t
    !******************************************************************************


















end module type_bcdata
