module type_element_info
    use mod_kinds,      only: ik
    use mod_constants,  only: NO_ID
    implicit none



    !> Container for locating a given element in the mesh data, via indices
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !---------------------------------------------------------------------
    type, public :: element_info_t

        integer(ik) :: idomain_g         = NO_ID
        integer(ik) :: idomain_l         = NO_ID
        integer(ik) :: ielement_g        = NO_ID
        integer(ik) :: ielement_l        = NO_ID
        integer(ik) :: iproc             = NO_ID
        integer(ik) :: pelem_ID          = NO_ID
        integer(ik) :: coordinate_system = NO_ID

        integer(ik) :: eqn_ID          = NO_ID
        integer(ik) :: nfields         = NO_ID
        integer(ik) :: ntime           = NO_ID
        integer(ik) :: nterms_s        = NO_ID
        integer(ik) :: nterms_c        = NO_ID
        integer(ik) :: dof_start       = NO_ID
        integer(ik) :: dof_local_start = NO_ID

        ! parallel access indices for native storage
        integer(ik) :: recv_comm    = NO_ID
        integer(ik) :: recv_domain  = NO_ID
        integer(ik) :: recv_element = NO_ID
        ! parallel access indices for petsc storage
        integer(ik) :: recv_dof     = NO_ID

    end type element_info_t
    !*********************************************************************

    interface element_info
        module procedure element_info_constructor
    end interface element_info



contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   5/7/2019
    !!
    !----------------------------------------------------------------------------------
    function element_info_constructor(idomain_g, idomain_l, ielement_g, ielement_l, iproc, pelem_ID, coordinate_system, eqn_ID, nfields, ntime, nterms_s, nterms_c, dof_start, dof_local_start, recv_comm, recv_domain, recv_element, recv_dof) result(elem_info)
        integer(ik),    intent(in)  :: idomain_g
        integer(ik),    intent(in)  :: idomain_l
        integer(ik),    intent(in)  :: ielement_g
        integer(ik),    intent(in)  :: ielement_l
        integer(ik),    intent(in)  :: iproc
        integer(ik),    intent(in)  :: pelem_ID
        integer(ik),    intent(in)  :: coordinate_system
        integer(ik),    intent(in)  :: eqn_ID
        integer(ik),    intent(in)  :: nfields
        integer(ik),    intent(in)  :: ntime
        integer(ik),    intent(in)  :: nterms_s
        integer(ik),    intent(in)  :: nterms_c
        integer(ik),    intent(in)  :: dof_start
        integer(ik),    intent(in)  :: dof_local_start
        integer(ik),    intent(in)  :: recv_comm
        integer(ik),    intent(in)  :: recv_domain
        integer(ik),    intent(in)  :: recv_element
        integer(ik),    intent(in)  :: recv_dof

        type(element_info_t)   :: elem_info

        elem_info%idomain_g         = idomain_g
        elem_info%idomain_l         = idomain_l
        elem_info%ielement_g        = ielement_g
        elem_info%ielement_l        = ielement_l
        elem_info%iproc             = iproc
        elem_info%pelem_ID          = pelem_ID
        elem_info%coordinate_system = coordinate_system
        elem_info%eqn_ID            = eqn_ID
        elem_info%nfields           = nfields
        elem_info%ntime             = ntime
        elem_info%nterms_s          = nterms_s
        elem_info%nterms_c          = nterms_c
        elem_info%dof_start         = dof_start
        elem_info%dof_local_start   = dof_local_start
        elem_info%recv_comm         = recv_comm
        elem_info%recv_domain       = recv_domain
        elem_info%recv_element      = recv_element
        elem_info%recv_dof          = recv_dof

    end function element_info_constructor
    !**********************************************************************************












end module type_element_info
