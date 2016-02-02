module type_chimera_donor



    !> A Chimera donor container for sending element data to Chimera receivers
    !! across domains.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-----------------------------------------------------------------------------------
    type, public :: chimera_donor_t

!        integer(ik), allocatable    :: ielem_send(:,:)      !< Element indices to be sent  (ielem_local, idomain_send)

    contains

!        procedure   :: send     !< Send element solution vectors to Chimera receivers

    end type chimera_donor_t
    !***********************************************************************************








end module type_chimera_donor
