module type_bc_group
    use type_bcvector,  only: bcvector_t
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/9/2016
    !!
    !-------------------------------------------------------------------------
    type, public :: bc_group_t
        
        character(:),       allocatable :: name         ! Boundary State Group name
        character(:),       allocatable :: family       ! Boundary State Group family
        type(bcvector_t)                :: bc_states    ! Vector of boundary condition state functions for each group.

    end type bc_group_t
    !*************************************************************************









end module type_bc_group
