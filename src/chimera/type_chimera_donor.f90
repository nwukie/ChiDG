module type_chimera_donor
    use mod_kinds,      only: ik
    use type_ivector,   only: ivector_t
    use mod_chidg_mpi,  only: IRANK
    implicit none



    !>  A Chimera donor container for sending element data to Chimera receivers
    !!  across domains.
    !!
    !!  Contains a list of elements being donated. On and Off-processor
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/12/2016
    !!
    !-----------------------------------------------------------------------------------
    type, public :: chimera_donor_t

        integer(ik) :: donor_domain_g
        integer(ik) :: donor_domain_l
        integer(ik) :: donor_element_g
        integer(ik) :: donor_element_l
        integer(ik) :: donor_proc

        ! Processors receiver faces are located on
        type(ivector_t) :: receiver_procs

    contains

        procedure   :: nrecipients

    end type chimera_donor_t
    !***********************************************************************************






contains




    !>
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2017
    !!
    !-----------------------------------------------------------------------------------
    function nrecipients(self) result(nrec)
        class(chimera_donor_t), intent(in)  :: self

        integer(ik) :: nrec

        nrec = self%receiver_procs%size()

    end function nrecipients
    !***********************************************************************************




!    !>  For a domain, add an element to the list of elements registered as Chimera donors.
!    !!
!    !!  @author Nathan A. Wukie (AFRL)
!    !!  @date   7/12/2016
!    !!
!    !!
!    !-------------------------------------------------------------------------------------------
!    subroutine add_donor(self,domain_g,domain_l,element_g,element_l,receiver_proc)
!        class(chimera_donor_t),     intent(inout)   :: self
!        integer(ik),                intent(in)      :: domain_g
!        integer(ik),                intent(in)      :: domain_l
!        integer(ik),                intent(in)      :: element_g
!        integer(ik),                intent(in)      :: element_l
!        integer(ik),                intent(in)      :: receiver_proc
!
!        integer(ik) :: idonor, ndonors
!        integer(ik) :: idonor_domain_g, idonor_element_g, ireceiver_proc
!        logical     :: already_added
!
!        ndonors = self%donor_domain_g%size()
!
!        !
!        ! Check if the donor was already added
!        !
!        already_added = .false.
!        do idonor = 1,ndonors
!
!            idonor_domain_g  = self%donor_domain_g%at(idonor)
!            idonor_element_g = self%donor_element_g%at(idonor)
!            ireceiver_proc   = self%receiver_proc%at(idonor)
!
!            already_added = ( (idonor_domain_g == domain_g)   .and. &
!                              (idonor_element_g == element_g) .and. &
!                              (ireceiver_proc == receiver_proc) )
!
!
!            if (already_added) then
!                exit
!            end if
!
!        end do
!
!
!
!        !
!        ! If we got all the way through the list above without exiting with an 'already_added' status, add the donor to the list.
!        !
!        if ( .not. already_added ) then
!            call self%donor_domain_g%push_back(domain_g)
!            call self%donor_domain_l%push_back(domain_l)
!            call self%donor_element_g%push_back(element_g)
!            call self%donor_element_l%push_back(element_l)
!            call self%receiver_proc%push_back(receiver_proc)
!        end if
!
!
!
!
!    end subroutine add_donor
!    !*******************************************************************************************
!
!
!
!
!
!
!    !>  Return the number of donor elements in the list
!    !!
!    !!  @author Nathan A. Wukie (AFRL)
!    !!  @date   7/13/2016
!    !!
!    !!
!    !!
!    !-------------------------------------------------------------------------------------------
!    function ndonors(self) result(res)
!        class(chimera_donor_t),     intent(in)  :: self
!
!        integer(ik) :: res
!
!        res = self%donor_domain_g%size()
!
!    end function ndonors
!    !********************************************************************************************
!
!
!
!
!
!
!    !>  Clear the list of registered Chimera donor elements.
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   3/11/2016
!    !!
!    !!
!    !-------------------------------------------------------------------------------------------
!    subroutine clear(self)
!        class(chimera_donor_t), intent(inout)   :: self
!
!        call self%donor_domain_g%clear()
!        call self%donor_domain_l%clear()
!        call self%donor_element_g%clear()
!        call self%donor_element_l%clear()
!        call self%receiver_proc%clear()
!
!    end subroutine clear
!    !*******************************************************************************************



end module type_chimera_donor
