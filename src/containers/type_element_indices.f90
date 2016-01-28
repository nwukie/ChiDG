module type_element_indices
    use mod_kinds,  only: ik



    !> Container for locating a given element in the mesh data, via indices
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !---------------------------------------------------------------------
    type, public :: element_indices_t

        integer(ik) :: idomain
        integer(ik) :: ielement

    end type element_indices_t
    !*********************************************************************




end module type_element_indices
