module type_bc_patch
#include <messenger.h>
    use mod_kinds,      only: ik
    use type_bc_coupled_element,    only: bc_coupled_element_t
    use type_ivector,   only: ivector_t
    implicit none


    !>  A boundary condition patch
    !!
    !!  This contains the boundary condition geometry description. A domain, 
    !!  element, and face index for each face the boundary condition is 
    !!  operating on.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!  @date   2/27/2017   ! updated for added generality
    !!
    !!
    !!
    !----------------------------------------------------------------------------
    type, public :: bc_patch_t

        integer(ik)                     :: patch_ID

        ! List of faces. Each combination(domain,element,face) in the
        ! vectors defines a face in the patch.
        type(ivector_t)                 :: idomain_g_
        type(ivector_t)                 :: idomain_l_
        type(ivector_t)                 :: ielement_g_
        type(ivector_t)                 :: ielement_l_
        type(ivector_t)                 :: iface_

        
        ! For each face in the patch, a list of elements it is coupled with
        type(bc_coupled_element_t), allocatable    :: coupling(:)
        !type(ivector_t),    allocatable :: idomain_g_coupled(:)
        !type(ivector_t),    allocatable :: idomain_l_coupled(:)
        !type(ivector_t),    allocatable :: ielement_g_coupled(:)
        !type(ivector_t),    allocatable :: ielement_l_coupled(:)
        !type(ivector_t),    allocatable :: proc_coupled(:)

        !type(ivector_t),    allocatable :: recv_comm(:)
        !type(ivector_t),    allocatable :: recv_domain(:)
        !type(ivector_t),    allocatable :: recv_element(:)

    contains

        ! Patch faces
        procedure   :: add_face
        procedure   :: nfaces

        ! Patch face coupling
        procedure   :: add_coupled_element
        procedure   :: ncoupled_elements
        procedure   :: set_coupled_element_recv

        ! Return indices for a given face_ID
        procedure   :: idomain_g
        procedure   :: idomain_l
        procedure   :: ielement_g
        procedure   :: ielement_l
        procedure   :: iface


        ! Parallal communication patterns
        procedure   :: get_recv_procs
        procedure   :: get_send_procs


    end type bc_patch_t
    !****************************************************************************






contains




    !>  Add a face to the boundary condition patch.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !----------------------------------------------------------------------------
    function add_face(self,idomain_g,idomain_l,ielement_g,ielement_l,iface) result(iface_bc)
        class(bc_patch_t),  intent(inout)   :: self
        integer(ik),        intent(in)      :: idomain_g
        integer(ik),        intent(in)      :: idomain_l
        integer(ik),        intent(in)      :: ielement_g
        integer(ik),        intent(in)      :: ielement_l
        integer(ik),        intent(in)      :: iface

        integer(ik)                     :: nfaces, ierr, iface_bc

        !
        ! Add face to bc_patch list
        !
        call self%idomain_g_%push_back(idomain_g)
        call self%idomain_l_%push_back(idomain_l)
        call self%ielement_g_%push_back(ielement_g)
        call self%ielement_l_%push_back(ielement_l)
        call self%iface_%push_back(iface)


        !
        ! Get location of face in bc_patch
        !
        iface_bc = self%iface_%size()


        !
        ! Extend coupling storage. Need a vector for each bc_face
        !
        if (allocated(self%idomain_g_coupled)) deallocate(self%idomain_g_coupled,   &
                                                          self%idomain_l_coupled,   &
                                                          self%ielement_g_coupled,  &
                                                          self%ielement_l_coupled,  &
                                                          self%proc_coupled,        &
                                                          self%recv_comm,           &
                                                          self%recv_domain,         &
                                                          self%recv_element)

        nfaces = self%nfaces()
        allocate(self%idomain_g_coupled(nfaces),    &
                 self%idomain_l_coupled(nfaces),    &
                 self%ielement_g_coupled(nfaces),   &
                 self%ielement_l_coupled(nfaces),   &
                 self%proc_coupled(nfaces),         &
                 self%recv_comm(nfaces),            &
                 self%recv_domain(nfaces),          &
                 self%recv_element(nfaces), stat=ierr)
        if (ierr /= 0) call AllocationError


    end function add_face
    !*****************************************************************************






    !>  Return the number of faces in the boundary condition patch
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !-----------------------------------------------------------------------------
    function nfaces(self) result(nfaces_bc)
        class(bc_patch_t),  intent(in)  :: self

        integer(ik) :: nfaces_bc

        nfaces_bc = self%iface_%size()

    end function nfaces
    !******************************************************************************






    !>  Return the global domain index of a boundary condition face
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2017
    !!
    !!
    !------------------------------------------------------------------------------
    function idomain_g(self,ind) result(idom_bc)
        class(bc_patch_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: ind

        integer(ik) :: idom_bc

        idom_bc = self%idomain_g_%at(ind)

    end function idomain_g
    !*******************************************************************************


    !>  Return the local domain index of a boundary condition face
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    function idomain_l(self,ind) result(idom_bc)
        class(bc_patch_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: ind

        integer(ik) :: idom_bc

        idom_bc = self%idomain_l_%at(ind)

    end function idomain_l
    !********************************************************************************




    !>  Return the global element index of a boundary condition face
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2017
    !!
    !!
    !--------------------------------------------------------------------------------
    function ielement_g(self,ind) result(ielem_bc)
        class(bc_patch_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: ind

        integer(ik) :: ielem_bc

        ielem_bc = self%ielement_g_%at(ind)

    end function ielement_g
    !*********************************************************************************



    !>  Return the local element index of a boundary condition face
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !---------------------------------------------------------------------------------
    function ielement_l(self,ind) result(ielem_bc)
        class(bc_patch_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: ind

        integer(ik) :: ielem_bc

        ielem_bc = self%ielement_l_%at(ind)

    end function ielement_l
    !**********************************************************************************







    !>  Return the face index of a boundary condition face
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !----------------------------------------------------------------------------------
    function iface(self,ind) result(iface_bc)
        class(bc_patch_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: ind

        integer(ik) :: iface_bc

        iface_bc = self%iface_%at(ind)

    end function iface
    !**********************************************************************************








    !>  For a face in the bc_patch, (iface_bc), add an element that it is coupled with.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2017
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine add_coupled_element(self,iface_bc,idomain_g, idomain_l, ielement_g, ielement_l, proc)
        class(bc_patch_t),  intent(inout)   :: self
        integer(ik),        intent(in)      :: iface_bc
        integer(ik),        intent(in)      :: idomain_g
        integer(ik),        intent(in)      :: idomain_l
        integer(ik),        intent(in)      :: ielement_g
        integer(ik),        intent(in)      :: ielement_l
        integer(ik),        intent(in)      :: proc

        integer(ik) :: idomain_g_coupled, ielement_g_coupled, ielem_coupled
        logical     :: already_added

        !
        ! Check if element has already been added to coupling list
        ! for the specified face
        !
        already_added = .false.
        do ielem_coupled = 1,self%ncoupled_elements(iface_bc)
            idomain_g_coupled  = self%idomain_g_coupled(iface_bc)%at(ielem_coupled)
            ielement_g_coupled = self%ielement_g_coupled(iface_bc)%at(ielem_coupled)

            if ( (idomain_g_coupled  == idomain_g ) .and. &
                 (ielement_g_coupled == ielement_g) ) already_added = .true.
            if (already_added) exit

        end do

        !
        ! If coupled element wasn't already added to the list, add all indices here
        !
        if (.not. already_added) then
            call self%idomain_g_coupled(iface_bc)%push_back(idomain_g)
            call self%idomain_l_coupled(iface_bc)%push_back(idomain_l)
            call self%ielement_g_coupled(iface_bc)%push_back(ielement_g)
            call self%ielement_l_coupled(iface_bc)%push_back(ielement_l)
            call self%proc_coupled(iface_bc)%push_back(proc)
        end if

    end subroutine add_coupled_element
    !***********************************************************************************







    !>  For a face in the bc_patch, (iface_bc), return the number of elements it is 
    !!  coupled with.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2017
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    function ncoupled_elements(self,iface_bc) result(res)
        class(bc_patch_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: iface_bc

        integer(ik) :: res

        res = self%ielement_g_coupled(iface_bc)%size()

    end function ncoupled_elements
    !************************************************************************************








    !>  Set the chidg_vector recv access data for a parallel coupled element.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/10/2017
    !!
    !!
    !----------------------------------------------------------------------------------
    subroutine set_coupled_element_recv(self,face_ID,idomain_g,ielement_g,recv_comm,recv_domain,recv_element)
        class(bc_patch_t),  intent(inout)   :: self
        integer(ik),        intent(in)      :: face_ID
        integer(ik),        intent(in)      :: idomain_g
        integer(ik),        intent(in)      :: ielement_g
        integer(ik),        intent(in)      :: recv_comm
        integer(ik),        intent(in)      :: recv_domain
        integer(ik),        intent(in)      :: recv_element

        character(:),   allocatable :: user_msg
        integer(ik)                 :: icoupled, loc
        logical                     :: element_found

        !
        ! Find the index associated with the element (idomain_g,ielement_g)
        !
        element_found = .false.
        do icoupled = 1,self%ncoupled_elements(face_ID)

            element_found = (idomain_g  = self%idomain_g_coupled(face_ID)%at(icoupled)) .and. &
                            (ielement_g = self%ielement_g_coupled(face_ID)%at(icoupled))
            if (element_found) loc = icoupled
            if (element_found) exit

        end do !icoupled

        user_msg = "bc_patch%set_coupled_element_recv: did not find element coupling."
        if (.not. element_found) call chidg_signal(FATAL,user_msg)


        self%recv_comm(face_ID)%data_(loc)    = recv_comm
        self%recv_domain(face_ID)%data_(loc)  = recv_domain
        self%recv_element(face_ID)%data_(loc) = recv_element

    end subroutine set_coupled_element_recv
    !***********************************************************************************






    !>  Return the processor ranks that the current bc_patch is receiving information
    !!  from.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/10/2017
    !!
    !!
    !------------------------------------------------------------------------------------
    function get_recv_procs(self) result(recv_procs_array)
        class(bc_patch_t),  intent(in)  :: self

        integer(ik),    allocatable :: recv_procs_array(:)
        type(ivector_t)             :: recv_procs
        integer(ik)                 :: face_ID, icoupled, proc
        logical                     :: comm_proc


        !
        ! Accumulate recv procs from each patch face
        !
        do face_ID = 1,self%nfaces()

            ! Loop through coupled elements, accumulate off-processor coupling ranks
            do icoupled = 1,self%ncoupled_elements(face_ID)
                proc = self%proc_coupled(face_ID)%at(icoupled)
                comm_proc = (proc /= IRANK)
                if (comm_proc) call recv_procs%push_back_unique(proc)
            end do !icoupled

        end do !face_ID


        !
        ! Return as integer array
        !
        recv_procs_array = recv_procs%data()

    end function get_recv_procs
    !*************************************************************************************







    !>  Return the processor ranks that the current bc_patch is sending information
    !!  to.
    !!
    !!  Currently, this pattern is the same as the receive pattern. So we just call
    !!  the routine: get_recv_procs
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/10/2017
    !!
    !!
    !------------------------------------------------------------------------------------
    function get_send_procs(self) result(send_procs_array)
        class(bc_patch_t),  intent(in)  :: self

        integer(ik),    allocatable :: send_procs_array(:)

        !
        ! Comm pattern currently same as recv. All processors we are receiving from, we 
        ! are also sending to.
        !
        send_procs_array = self%get_recv_procs()

    end function get_send_procs
    !*************************************************************************************









end module type_bc_patch
