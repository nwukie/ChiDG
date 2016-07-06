module type_seed
    use mod_kinds,  only: ik
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

    end type seed_t
    !**********************************************************************************************





end module type_seed
