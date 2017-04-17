module type_bc_patch
#include <messenger.h>
    use mod_kinds,                  only: ik
    use type_ivector,               only: ivector_t
    use type_bc_element_coupling,   only: bc_element_coupling_t
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
        type(bc_element_coupling_t), allocatable    :: coupling(:)

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


        ! Parallel communication patterns
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

        integer(ik) :: nfaces, ierr, iface_bc

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
        if (allocated(self%coupling)) deallocate(self%coupling)

        nfaces = self%nfaces()
        allocate(self%coupling(nfaces), stat=ierr)
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
    function idomain_g(self,face_ID) result(idomain_g_)
        class(bc_patch_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: face_ID

        integer(ik) :: idomain_g_

        idomain_g_ = self%idomain_g_%at(face_ID)

    end function idomain_g
    !*******************************************************************************


    !>  Return the local domain index of a boundary condition face
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    function idomain_l(self,face_ID) result(idomain_l_)
        class(bc_patch_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: face_ID

        integer(ik) :: idomain_l_

        idomain_l_ = self%idomain_l_%at(face_ID)

    end function idomain_l
    !********************************************************************************




    !>  Return the global element index of a boundary condition face
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2017
    !!
    !!
    !--------------------------------------------------------------------------------
    function ielement_g(self,face_ID) result(ielement_g_)
        class(bc_patch_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: face_ID

        integer(ik) :: ielement_g_

        ielement_g_ = self%ielement_g_%at(face_ID)

    end function ielement_g
    !*********************************************************************************



    !>  Return the local element index of a boundary condition face
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !---------------------------------------------------------------------------------
    function ielement_l(self,face_ID) result(ielement_l_)
        class(bc_patch_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: face_ID

        integer(ik) :: ielement_l_

        ielement_l_ = self%ielement_l_%at(face_ID)

    end function ielement_l
    !**********************************************************************************







    !>  Return the face index of a boundary condition face
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !----------------------------------------------------------------------------------
    function iface(self,face_ID) result(iface_)
        class(bc_patch_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: face_ID

        integer(ik) :: iface_

        iface_ = self%iface_%at(face_ID)

    end function iface
    !**********************************************************************************








    !>  For a face in the bc_patch, (face_ID), add an element that it is coupled with.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2017
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine add_coupled_element(self,face_ID,idomain_g, idomain_l, ielement_g, ielement_l, iface, proc)
        class(bc_patch_t),  intent(inout)   :: self
        integer(ik),        intent(in)      :: face_ID
        integer(ik),        intent(in)      :: idomain_g
        integer(ik),        intent(in)      :: idomain_l
        integer(ik),        intent(in)      :: ielement_g
        integer(ik),        intent(in)      :: ielement_l
        integer(ik),        intent(in)      :: iface
        integer(ik),        intent(in)      :: proc


        call self%coupling(face_ID)%add_coupled_element(idomain_g,idomain_l,ielement_g,ielement_l,iface,proc)


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

        res = self%coupling(iface_bc)%ncoupled_elements()

    end function ncoupled_elements
    !************************************************************************************








    !>  Set the chidg_vector recv access data for a parallel coupled element.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/10/2017
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine set_coupled_element_recv(self,face_ID,idomain_g,ielement_g,recv_comm,recv_domain,recv_element)
        class(bc_patch_t),  intent(inout)   :: self
        integer(ik),        intent(in)      :: face_ID
        integer(ik),        intent(in)      :: idomain_g
        integer(ik),        intent(in)      :: ielement_g
        integer(ik),        intent(in)      :: recv_comm
        integer(ik),        intent(in)      :: recv_domain
        integer(ik),        intent(in)      :: recv_element


        call self%coupling(face_ID)%set_coupled_element_recv(idomain_g,ielement_g,recv_comm,recv_domain,recv_element)


    end subroutine set_coupled_element_recv
    !************************************************************************************






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
                proc = self%coupling(face_ID)%proc(icoupled)
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
