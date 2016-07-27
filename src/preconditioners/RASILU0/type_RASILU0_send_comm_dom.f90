module type_RASILU0_send_comm_dom
#include <messenger.h>
    use mod_kinds,      only: ik
    use type_ivector,   only: ivector_t
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/25/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------
    type, public :: RASILU0_send_comm_dom_t

        integer(ik)                     :: idomain_g
        integer(ik)                     :: idomain_l
        type(ivector_t)                 :: elem_send    !< For the current domain, list of elements to send
        type(ivector_t),    allocatable :: blk_send(:)  !< For each element, a list of block indices to send

    contains

!        procedure   :: init

    end type RASILU0_send_comm_dom_t
    !**************************************************************************************





contains





!    !>
!    !!
!    !!  @author Nathan A. Wukie (AFRL)
!    !!  @date   7/25/2016
!    !!
!    !!
!    !!
!    !----------------------------------------------------------------------------------------
!    subroutine init(self,idomain_g,idomain_l,nelem_send)
!        class(RASILU0_send_comm_dom_t), intent(inout)   :: self
!        integer(ik),                    intent(in)      :: idomain_g
!        integer(ik),                    intent(in)      :: idomain_l
!        integer(ik),                    intent(in)      :: nelem_send
!
!        integer(ik) :: ierr
!
!
!        self%idomain_g = idomain_g
!        self%idomain_l = idomain_l
!
!        allocate(self%blk_send(nelem_send), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!    end subroutine init
!    !*****************************************************************************************










end module type_RASILU0_send_comm_dom
