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

        integer(ik) :: idomain_g
        integer(ik) :: idomain_l
        integer(ik) :: ielement_g
        integer(ik) :: ielement_l
        integer(ik) :: iproc

        integer(ik) :: neqns
        integer(ik) :: nterms_s
        integer(ik) :: nterms_c

    end type element_indices_t
    !*********************************************************************




end module type_element_indices
