module type_bc_coupling
    use mod_kinds,      only: ik
    use type_ivector,   only: ivector_t
    implicit none



    !>  Container for holding coupling information for a boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/16/2016
    !!
    !!
    !-----------------------------------------------------------------------------
    type, public :: bc_coupling_t

        integer(ik),        allocatable :: elems(:)             !< Array of element indices associated with the boundary condition.
        type(ivector_t),    allocatable :: coupled_elems(:)     !< For each element in self%elems(:), a vector of coupled elements.

    end type bc_coupling_t
    !*****************************************************************************





end module type_bc_coupling
