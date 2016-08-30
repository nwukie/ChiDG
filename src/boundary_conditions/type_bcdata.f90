module type_bcdata
    use mod_kinds,                  only: ik
!    use type_bcwrapper,             only: bcwrapper_t
    use type_bc,                    only: bc_t
    use type_boundary_connectivity, only: boundary_connectivity_t
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

        character(len=:),               allocatable :: domain_              !< Domain name the bcdata is associated with
        !type(bcwrapper_t),              allocatable :: bcs(:)               !< Array of boundary conditions for each face
        type(bc_t),                     allocatable :: bcs(:)               !< Array of boundary conditions for each face
        type(boundary_connectivity_t),  allocatable :: bc_connectivity(:)   !< 
        !type(connectivity_t),   allocatable :: bc_connectivity(:)          !<

    contains


    end type bcdata_t
    !******************************************************************************


















end module type_bcdata
