module type_seed
    use mod_kinds,  only: ik



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
    !----------------------------------------------------------------------
    type, public :: seed_t

        integer(ik) :: idom
        integer(ik) :: ielem

    end type seed_t
    !**********************************************************************








end module type_seed
