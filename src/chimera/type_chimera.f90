module type_chimera
#include <messenger.h>
    use mod_kinds,              only: ik, rk
    use mod_constants,          only: NO_ID
    use type_chimera_receiver,  only: chimera_receiver_t, chimera_receiver
    use type_chimera_send,      only: chimera_send_t
    implicit none



    !> Main interface and container for Chimera data and operations.
    !! Holds chimera send/receive sets which are used to facilitate inter-domain communication
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !---------------------------------------------------------------------------------------------
    type, public :: chimera_t
    
        type(chimera_receiver_t),   allocatable :: recv(:)  ! receiver face instances. One for each Chimera receiver on a domain
        type(chimera_send_t),       allocatable :: send(:)  ! one chimera send instance for each donor on a domain being sent off-processor

    contains


        procedure   :: add_receiver
        procedure   :: new_receiver
        procedure   :: nreceivers
        procedure   :: find_receiver


        procedure   :: new_send
        procedure   :: find_send
        procedure   :: nsend

        procedure   :: clear

    end type chimera_t
    !*********************************************************************************************


contains







    !>  Add chimera receiver face to the chimera collection.
    !!
    !!  If receiver face matching the incoming face already exists, return its ID instead of
    !!  creating a new instance.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2016
    !!
    !----------------------------------------------------------------------------------------
    function add_receiver(self, idomain_g, idomain_l, ielement_g, ielement_l, iface, nnodes) result(ChiID)
        class(chimera_t),   intent(inout)   :: self
        integer(ik),        intent(in)      :: idomain_g
        integer(ik),        intent(in)      :: idomain_l
        integer(ik),        intent(in)      :: ielement_g
        integer(ik),        intent(in)      :: ielement_l
        integer(ik),        intent(in)      :: iface
        integer(ik),        intent(in)      :: nnodes

        integer(ik) :: ChiID


        !
        ! Check if receiver matching the incoming face already exists
        !
        if ( self%find_receiver(idomain_g, ielement_g, iface) == NO_ID ) then
            ChiID = self%new_receiver()
            self%recv(ChiID) = chimera_receiver(idomain_g, idomain_l, ielement_g, ielement_l, iface, nnodes)

        else
            ChiID = self%find_receiver(idomain_g,ielement_g,iface)
        end if

    end function add_receiver
    !*****************************************************************************************







    !>  Extend allocation of receiver instances, return identifier for new instance.
    !!
    !!  Returns recv_ID for chimera_receiver_t in domain%chimera%recv(recv_ID)
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2017
    !!
    !-----------------------------------------------------------------------------------------
    function new_receiver(self) result(recv_ID)
        class(chimera_t),   intent(inout)   :: self


        type(chimera_receiver_t),   allocatable :: temp(:)
        integer(ik)                             :: ierr, recv_ID



        ! Resize array storage
        allocate(temp(self%nreceivers() + 1), stat=ierr)


        ! Copy previously initialized instances to new array.
        if (self%nreceivers() > 0) then
            temp(1:size(self%recv)) = self%recv(1:size(self%recv))
        end if


        ! Move resized temp allocation back to parent. 
        call move_alloc(temp,self%recv)


        ! Set domain identifier of newly allocated domain that will be returned
        recv_ID = self%nreceivers()


    end function new_receiver
    !*****************************************************************************************







    !>  Return number of receiver faces currently exist in the chimera collection.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2017
    !!
    !-----------------------------------------------------------------------------------------
    function nreceivers(self) result(nrecv)
        class(chimera_t),   intent(in)  :: self

        integer(ik) :: nrecv

        if (allocated(self%recv)) then
            nrecv = size(self%recv)
        else
            nrecv = 0
        end if

    end function nreceivers
    !*****************************************************************************************




    !>  Return number of send elements currently exist in the chimera collection.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/25/2017
    !!
    !-----------------------------------------------------------------------------------------
    function nsend(self) result(nsend_)
        class(chimera_t),   intent(in)  :: self

        integer(ik) :: nsend_

        if (allocated(self%send)) then
            nsend_ = size(self%send)
        else
            nsend_ = 0
        end if

    end function nsend
    !*****************************************************************************************






    !>  Given a global face location, return the index of a chimera_receiver instance 
    !!  corresponding to the face, if one exists in the chimera collection.
    !!
    !!  If no chimera_receiver instance can be found that matches the incoming description
    !!  the return NO_ID.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2017
    !!
    !------------------------------------------------------------------------------------------
    function find_receiver(self,idomain_g,ielement_g,iface) result(recv_ID)
        class(chimera_t),   intent(in)  :: self
        integer(ik),        intent(in)  :: idomain_g
        integer(ik),        intent(in)  :: ielement_g
        integer(ik),        intent(in)  :: iface

        integer(ik) :: recv_ID, irecv 

        recv_ID = NO_ID
        do irecv = 1,self%nreceivers()
            if ( (self%recv(irecv)%idomain_g  == idomain_g ) .and. &
                 (self%recv(irecv)%ielement_g == ielement_g) .and. &
                 (self%recv(irecv)%iface      == iface     ) ) then
                 recv_ID = irecv
                 exit
            end if
        end do !irecv

    end function find_receiver
    !******************************************************************************************





    !>  Given a global element location, return the index of a chimera_send instance 
    !!  corresponding to the element, if one exists in the chimera collection.
    !!
    !!  If no chimera_send instance can be found that matches the incoming description
    !!  the return NO_ID.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/25/2017
    !!
    !------------------------------------------------------------------------------------------
    function find_send(self,idomain_g,ielement_g) result(send_ID)
        class(chimera_t),   intent(in)  :: self
        integer(ik),        intent(in)  :: idomain_g
        integer(ik),        intent(in)  :: ielement_g

        integer(ik) :: send_ID, isend

        send_ID = NO_ID
        do isend = 1,self%nsend()
            if ( (self%send(isend)%idomain_g == idomain_g) .and. &
                 (self%send(isend)%ielement_g == ielement_g) ) then
                 send_ID = isend
                 exit
            end if
        end do !isend

    end function find_send
    !******************************************************************************************







    !>  Extend allocation of send instances, return identifier for new instance.
    !!
    !!  Returns send_ID for chimera_send_t in domain%chimera%send(send_ID)
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/25/2017
    !!
    !-----------------------------------------------------------------------------------------
    function new_send(self) result(send_ID)
        class(chimera_t),   intent(inout)   :: self


        type(chimera_send_t),  allocatable :: temp(:)
        integer(ik)                         :: ierr, send_ID



        ! Resize array storage
        allocate(temp(self%nsend() + 1), stat=ierr)


        ! Copy previously initialized instances to new array.
        if (self%nsend() > 0) then
            temp(1:size(self%send)) = self%send(1:size(self%send))
        end if


        ! Move resized temp allocation back to parent. 
        call move_alloc(temp,self%send)


        ! Set domain identifier of newly allocated domain that will be returned
        send_ID = self%nsend()


    end function new_send
    !*****************************************************************************************






    !>
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(chimera_t),   intent(inout)  :: self

        if (allocated(self%recv )) deallocate(self%recv )
        if (allocated(self%send )) deallocate(self%send )

    end subroutine clear
    !******************************************************************************************



end module type_chimera
