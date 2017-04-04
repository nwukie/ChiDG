module type_element_info
    use mod_kinds,  only: ik



    !> Container for locating a given element in the mesh data, via indices
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !---------------------------------------------------------------------
    type, public :: element_info_t

        integer(ik) :: idomain_g
        integer(ik) :: idomain_l
        integer(ik) :: ielement_g
        integer(ik) :: ielement_l
        integer(ik) :: iproc

        integer(ik) :: eqn_ID
        integer(ik) :: neqns
        integer(ik) :: nterms_s
        integer(ik) :: nterms_c

    end type element_info_t
    !*********************************************************************




end module type_element_info
