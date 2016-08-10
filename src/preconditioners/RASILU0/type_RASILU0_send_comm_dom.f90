module type_RASILU0_send_comm_dom
#include <messenger.h>
    use mod_kinds,      only: ik
    use type_ivector,   only: ivector_t
    implicit none




    !>  Container for storing the overlap element information for a domain
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

    end type RASILU0_send_comm_dom_t
    !**************************************************************************************





contains







end module type_RASILU0_send_comm_dom
