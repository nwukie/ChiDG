module type_seed
    use mod_kinds,      only: ik
    use mod_constants,  only: NO_PROC, NO_ID
    implicit none


    !> Container that holds information on the element solution being 
    !!  linearized with respect to.
    !!
    !!  For example, if we were computing:
    !!
    !!  \f$     \frac{\partial F}{\partial Q_{idom,ielem}}      \f$
    !!
    !!  This container stores the indices of idom,ielem so the correct
    !!  solution variables are initialized in the automatic differentiation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !-------------------------------------------------------------------------------------
    type, public :: seed_t

        integer(ik) :: idomain_g
        integer(ik) :: idomain_l
        integer(ik) :: ielement_g
        integer(ik) :: ielement_l
        integer(ik) :: nfields
        integer(ik) :: nterms_s
        integer(ik) :: nnodes_r
        integer(ik) :: iproc
        integer(ik) :: itime
        integer(ik) :: dof_start
        integer(ik) :: xdof_start

        ! If seed is on another processor, these are its location in the recv 
        ! container on the current processor. Otherwise, not used.
        integer(ik) :: recv_comm
        integer(ik) :: recv_domain
        integer(ik) :: recv_element

    contains

        procedure   :: init
        procedure   :: clear

    end type seed_t
    !*************************************************************************************






contains



    !>  Initialize the seed data components.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/6/2016
    !!
    !-------------------------------------------------------------------------------------
    subroutine init(self,idomain_g,idomain_l,ielement_g,ielement_l,nfields,nterms_s,nnodes_r,iproc,itime,dof_start,xdof_start,recv_comm,recv_domain,recv_element)
        class(seed_t),  intent(inout)   :: self
        integer(ik),    intent(in)      :: idomain_g
        integer(ik),    intent(in)      :: idomain_l
        integer(ik),    intent(in)      :: ielement_g
        integer(ik),    intent(in)      :: ielement_l
        integer(ik),    intent(in)      :: nfields
        integer(ik),    intent(in)      :: nterms_s
        integer(ik),    intent(in)      :: nnodes_r
        integer(ik),    intent(in)      :: iproc
        integer(ik),    intent(in)      :: itime
        integer(ik),    intent(in)      :: dof_start
        integer(ik),    intent(in)      :: xdof_start
        integer(ik),    intent(in)      :: recv_comm
        integer(ik),    intent(in)      :: recv_domain
        integer(ik),    intent(in)      :: recv_element

        self%idomain_g    = idomain_g
        self%idomain_l    = idomain_l
        self%ielement_g   = ielement_g
        self%ielement_l   = ielement_l
        self%nfields      = nfields
        self%nterms_s     = nterms_s
        self%nnodes_r     = nnodes_r
        self%iproc        = iproc
        self%itime        = itime
        self%dof_start    = dof_start
        self%xdof_start   = xdof_start

        self%recv_comm    = recv_comm
        self%recv_domain  = recv_domain
        self%recv_element = recv_element

    end subroutine init
    !**************************************************************************************




    !>  Clear the seed data components.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/13/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine clear(self)
        class(seed_t),  intent(inout)   :: self

        self%idomain_g    = 0
        self%idomain_l    = 0
        self%ielement_g   = 0
        self%ielement_l   = 0
        self%nfields      = 0
        self%nterms_s     = 0
        self%nnodes_r     = 0
        self%iproc        = NO_PROC
        self%itime        = 0
        self%dof_start    = NO_ID
        self%xdof_start   = NO_ID

        self%recv_comm    = 0
        self%recv_domain  = 0
        self%recv_element = 0

    end subroutine clear
    !**************************************************************************************











end module type_seed
