module type_chimera
#include <messenger.h>
    use mod_kinds,              only: ik, rk
    use mod_constants,          only: NO_ID
    use type_chimera_receiver,  only: chimera_receiver_t
    use type_chimera_donor,     only: chimera_donor_t
    implicit none



    !> Main interface and container for Chimera data and operations.
    !! Holds chimera send/receive sets which are used to facilitate inter-domain communication
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !---------------------------------------------------------------------------------------------
    type, public :: chimera_t
    
        !type(chimera_receiver_t)    :: recv
        !type(chimera_donor_t)       :: send

        type(chimera_receiver_t),   allocatable :: recv(:)
        type(chimera_donor_t),      allocatable :: donor(:)


    contains


        procedure   :: add_receiver
        procedure   :: new_receiver
        procedure   :: nreceivers
        procedure   :: find_receiver

        procedure   :: add_donor
        procedure   :: new_donor
        procedure   :: ndonors
        procedure   :: find_donor

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
    function add_receiver(self, idomain_g, idomain_l, ielement_g, ielement_l, iface, iproc) result(recv_ID)
        class(chimera_t),   intent(inout)   :: self
        integer(ik),        intent(in)      :: idomain_g
        integer(ik),        intent(in)      :: idomain_l
        integer(ik),        intent(in)      :: ielement_g
        integer(ik),        intent(in)      :: ielement_l
        integer(ik),        intent(in)      :: iface
        integer(ik),        intent(in)      :: iproc

        integer(ik) :: recv_ID


        !
        ! Check if receiver matching the incoming face already exists
        !
        if ( self%find_receiver(idomain_g, ielement_g, iface) == NO_ID ) then
            recv_ID = self%new_receiver()

            self%recv(recv_ID)%receiver_domain_g  = idomain_g
            self%recv(recv_ID)%receiver_domain_l  = idomain_l
            self%recv(recv_ID)%receiver_element_g = ielement_g
            self%recv(recv_ID)%receiver_element_l = ielement_l
            self%recv(recv_ID)%receiver_face      = iface
            self%recv(recv_ID)%receiver_proc      = iproc

        else
            recv_ID = self%find_receiver(idomain_g,ielement_g,iface)
        end if

    end function add_receiver
    !*****************************************************************************************




    !>  Add chimera donor element to the chimera collection.
    !!
    !!  If donor element matching the incoming element already exists, return its ID instead of
    !!  creating a new instance.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2016
    !!
    !----------------------------------------------------------------------------------------
    function add_donor(self, donor_domain_g, donor_domain_l, donor_element_g, donor_element_l, donor_proc) result(donor_ID)
        class(chimera_t),   intent(inout)   :: self
        integer(ik),        intent(in)      :: donor_domain_g
        integer(ik),        intent(in)      :: donor_domain_l
        integer(ik),        intent(in)      :: donor_element_g
        integer(ik),        intent(in)      :: donor_element_l
        integer(ik),        intent(in)      :: donor_proc

        integer(ik) :: donor_ID


        !
        ! Check if receiver matching the incoming face already exists
        !
        if ( self%find_donor(donor_domain_g, donor_element_g) == NO_ID ) then
            donor_ID = self%new_donor()

            self%donor(donor_ID)%donor_domain_g  = donor_domain_g
            self%donor(donor_ID)%donor_domain_l  = donor_domain_l
            self%donor(donor_ID)%donor_element_g = donor_element_g
            self%donor(donor_ID)%donor_element_l = donor_element_l
            self%donor(donor_ID)%donor_proc      = donor_proc

        else
            donor_ID = self%find_donor(donor_domain_g,donor_element_g)
        end if

    end function add_donor
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



    !>  Extend allocation of donor instances, return identifier for new instance.
    !!
    !!  Returns donor_ID for chimera_donor_t in domain%chimera%donor(donor_ID)
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2017
    !!
    !-----------------------------------------------------------------------------------------
    function new_donor(self) result(donor_ID)
        class(chimera_t),   intent(inout)   :: self


        type(chimera_donor_t),  allocatable :: temp(:)
        integer(ik)                         :: ierr, donor_ID



        ! Resize array storage
        allocate(temp(self%ndonors() + 1), stat=ierr)


        ! Copy previously initialized instances to new array.
        if (self%ndonors() > 0) then
            temp(1:size(self%donor)) = self%donor(1:size(self%donor))
        end if


        ! Move resized temp allocation back to parent. 
        call move_alloc(temp,self%donor)


        ! Set domain identifier of newly allocated domain that will be returned
        donor_ID = self%ndonors()


    end function new_donor
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



    !>  Return number of donor elements currently exist in the chimera collection.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2017
    !!
    !-----------------------------------------------------------------------------------------
    function ndonors(self) result(ndonors_)
        class(chimera_t),   intent(in)  :: self

        integer(ik) :: ndonors_

        if (allocated(self%donor)) then
            ndonors_ = size(self%donor)
        else
            ndonors_ = 0
        end if

    end function ndonors
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
            if ( (self%recv(irecv)%receiver_domain_g == idomain_g) .and. &
                 (self%recv(irecv)%receiver_element_g == ielement_g) .and. &
                 (self%recv(irecv)%receiver_face == iface) ) then
                 recv_ID = irecv
            end if
        end do !irecv

    end function find_receiver
    !******************************************************************************************








    !>  Given a global element location, return the index of a chimera_donor instance 
    !!  corresponding to the element, if one exists in the chimera collection.
    !!
    !!  If no chimera_donor instance can be found that matches the incoming description
    !!  the return NO_ID.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2017
    !!
    !------------------------------------------------------------------------------------------
    function find_donor(self,idomain_g,ielement_g) result(donor_ID)
        class(chimera_t),   intent(in)  :: self
        integer(ik),        intent(in)  :: idomain_g
        integer(ik),        intent(in)  :: ielement_g

        integer(ik) :: donor_ID, idonor

        donor_ID = NO_ID
        do idonor = 1,self%ndonors()
            if ( (self%donor(idonor)%donor_domain_g == idomain_g) .and. &
                 (self%donor(idonor)%donor_element_g == ielement_g) ) then
                 donor_ID = idonor
            end if
        end do !idonor

    end function find_donor
    !******************************************************************************************




    !>
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(chimera_t),   intent(inout)  :: self

        if (allocated(self%recv)) deallocate(self%recv)
        if (allocated(self%donor)) deallocate(self%donor)

    end subroutine clear
    !******************************************************************************************



end module type_chimera
