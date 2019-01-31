module type_element_info
    use mod_kinds,      only: ik
    use mod_constants,  only: NO_ID



    !> Container for locating a given element in the mesh data, via indices
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !---------------------------------------------------------------------
    type, public :: element_info_t

        integer(ik) :: idomain_g  = NO_ID
        integer(ik) :: idomain_l  = NO_ID
        integer(ik) :: ielement_g = NO_ID
        integer(ik) :: ielement_l = NO_ID
        integer(ik) :: iproc      = NO_ID
        integer(ik) :: pelem_ID   = NO_ID

        integer(ik) :: eqn_ID     = NO_ID
        integer(ik) :: nfields    = NO_ID
        integer(ik) :: nterms_s   = NO_ID
        integer(ik) :: nterms_c   = NO_ID
        integer(ik) :: dof_start  = NO_ID

        ! access indices in chidg_vector
        integer(ik) :: recv_comm    = NO_ID
        integer(ik) :: recv_domain  = NO_ID
        integer(ik) :: recv_element = NO_ID

    end type element_info_t
    !*********************************************************************




end module type_element_info
