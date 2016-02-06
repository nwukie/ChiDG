module type_bcdata
    use mod_kinds,          only: ik
    use type_bcwrapper,     only: bcwrapper_t
    implicit none





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------
    type, public :: bcdata_t

        character(len=:),   allocatable :: domain_      !< Domain name the bcdata is associated with
        type(bcwrapper_t),  allocatable :: bcs(:)       !< Array of boundary conditions for each face
        integer(ik),        allocatable :: bcface(:)    !< bc face index. XI_MIN, XI_MAX, ETA_MIN, etc.

    contains


    end type bcdata_t
    !******************************************************************************


















end module type_bcdata
