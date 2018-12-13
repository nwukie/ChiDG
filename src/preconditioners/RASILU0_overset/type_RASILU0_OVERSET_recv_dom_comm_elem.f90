module type_RASILU0_OVERSET_recv_dom_comm_elem
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use type_densematrix,   only: densematrix_t
    use type_ivector,       only: ivector_t
    implicit none



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    type, public :: RASILU0_OVERSET_recv_dom_comm_elem_t

        type(densematrix_t),    allocatable :: blks(:)          !< For a given overlapping element, the blocks being received.
        integer(ik),            allocatable :: blk_indices(:)   !< Indices of the blocks being received on their host proc
        integer(ik),            allocatable :: trans_elem(:)    !< Element index of block in transposed location
        integer(ik),            allocatable :: trans_blk(:)     !< Block index of block in transposed location

        
        type(ivector_t) :: lower                                !< indices of lower blocks
        type(ivector_t) :: upper                                !< indices of upper blocks
        type(ivector_t) :: diag                                 !< index of diagonal block. Should be size 1


    contains

        procedure   :: init

    end type RASILU0_OVERSET_recv_dom_comm_elem_t
    !***********************************************************************************************




contains




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine init(self,nblocks)
        class(RASILU0_OVERSET_recv_dom_comm_elem_t),    intent(inout)   :: self
        integer(ik),                            intent(in)      :: nblocks

        integer(ik) :: ierr

        allocate(self%blks(nblocks),            &
                 self%blk_indices(nblocks),     &
                 self%trans_elem(nblocks),      &
                 self%trans_blk(nblocks), stat=ierr)
        if (ierr /= 0) call AllocationError


    end subroutine init
    !***********************************************************************************************



end module type_RASILU0_OVERSET_recv_dom_comm_elem
