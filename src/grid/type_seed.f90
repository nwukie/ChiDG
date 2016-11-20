module type_seed
    use mod_kinds,      only: ik
    use mod_constants,  only: NO_PROC
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
    !----------------------------------------------------------------------------------------------
    type, public :: seed_t

        integer(ik) :: idomain_g
        integer(ik) :: idomain_l
        integer(ik) :: ielement_g
        integer(ik) :: ielement_l
        integer(ik) :: iproc

        ! If seed is on another processor, these are its location in the recv container on the current processor
        ! Otherwise, not used.
        integer(ik) :: recv_comm
        integer(ik) :: recv_domain
        integer(ik) :: recv_element

    contains

        procedure   :: init

        procedure   :: clear

    end type seed_t
    !**********************************************************************************************






contains



    !>  Initialize seed_t container.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/20/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init(self,idomain_g,idomain_l,ielement_g,ielement_l,iproc)
        class(seed_t),  intent(inout)   :: self
        integer(ik),    intent(in)      :: idomain_g
        integer(ik),    intent(in)      :: idomain_l
        integer(ik),    intent(in)      :: ielement_g
        integer(ik),    intent(in)      :: ielement_l
        integer(ik),    intent(in)      :: iproc

        self%idomain_g  = idomain_g
        self%idomain_l  = idomain_l
        self%ielement_g = ielement_g
        self%ielement_l = ielement_l
        self%iproc      = iproc

    end subroutine init
    !*********************************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/13/2016
    !!
    !!
    !----------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(seed_t),  intent(inout)   :: self

        self%idomain_g    = 0
        self%idomain_l    = 0
        self%ielement_g   = 0
        self%ielement_l   = 0
        self%iproc        = NO_PROC

        self%recv_comm    = 0
        self%recv_domain  = 0
        self%recv_element = 0

    end subroutine clear
    !**********************************************************************************************











end module type_seed
