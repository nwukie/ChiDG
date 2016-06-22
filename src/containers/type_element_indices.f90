module type_element_indices
    use mod_kinds,  only: ik



    !> Container for locating a given element in the mesh data, via indices
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !---------------------------------------------------------------------
    type, public :: element_indices_t

!        integer(ik) :: idomain
!        integer(ik) :: ielement
        integer(ik) :: idomain_g
        integer(ik) :: idomain_l
        integer(ik) :: ielement_g
        integer(ik) :: ielement_l

    end type element_indices_t
    !*********************************************************************




end module type_element_indices
