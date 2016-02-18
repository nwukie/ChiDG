module type_bcset_coupling
    use type_bc_coupling,   only: bc_coupling_t
    implicit none
        



    !>  Holds an array of bc_coupling_t instances describing the coupling of a 
    !!  boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/16/2016
    !!
    !!
    !----------------------------------------------------------------------------
    type, public :: bcset_coupling_t

        type(bc_coupling_t),    allocatable :: bc(:)

    end type bcset_coupling_t
    !****************************************************************************





end module type_bcset_coupling
