module type_bc_patch
#include <messenger.h>
    use mod_kinds,      only: ik
    use type_ivector,   only: ivector_t
    implicit none


    !>  A boundary condition patch
    !!
    !!  This contains the boundary condition geometry description. A domain, element, and face
    !!  index for each face the boundary condition is operating on.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!  @date   2/27/2017   ! updated for added generality
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
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
        type(ivector_t),    allocatable :: idomain_g_coupled(:)
        type(ivector_t),    allocatable :: idomain_l_coupled(:)
        type(ivector_t),    allocatable :: ielement_g_coupled(:)
        type(ivector_t),    allocatable :: ielement_l_coupled(:)
        type(ivector_t),    allocatable :: proc_coupled(:)

    contains

        procedure   :: add_face
        procedure   :: nfaces

        ! Return indices for a given bc face
        procedure   :: idomain_g
        procedure   :: idomain_l
        procedure   :: ielement_g
        procedure   :: ielement_l
        procedure   :: iface

        procedure   :: add_coupled_element
        procedure   :: ncoupled_elements

    end type bc_patch_t
    !******************************************************************************************






contains




    !>  Add a face to the boundary condition patch.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
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
                                                          self%proc_coupled)

        nfaces = self%nfaces()
        allocate(self%idomain_g_coupled(nfaces),    &
                 self%idomain_l_coupled(nfaces),    &
                 self%ielement_g_coupled(nfaces),   &
                 self%ielement_l_coupled(nfaces),   &
                 self%proc_coupled(nfaces), stat=ierr)
        if (ierr /= 0) call AllocationError


    end function add_face
    !******************************************************************************************






    !>  Return the number of faces in the boundary condition patch
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !------------------------------------------------------------------------------------------
    function nfaces(self) result(nfaces_bc)
        class(bc_patch_t),  intent(in)  :: self

        integer(ik) :: nfaces_bc

        nfaces_bc = self%iface_%size()

    end function nfaces
    !*******************************************************************************************






    !>  Return the global domain index of a boundary condition face
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    function idomain_g(self,ind) result(idom_bc)
        class(bc_patch_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: ind

        integer(ik) :: idom_bc

        idom_bc = self%idomain_g_%at(ind)

    end function idomain_g
    !******************************************************************************************


    !>  Return the local domain index of a boundary condition face
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    function idomain_l(self,ind) result(idom_bc)
        class(bc_patch_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: ind

        integer(ik) :: idom_bc

        idom_bc = self%idomain_l_%at(ind)

    end function idomain_l
    !******************************************************************************************




    !>  Return the global element index of a boundary condition face
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    function ielement_g(self,ind) result(ielem_bc)
        class(bc_patch_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: ind

        integer(ik) :: ielem_bc

        ielem_bc = self%ielement_g_%at(ind)

    end function ielement_g
    !******************************************************************************************



    !>  Return the local element index of a boundary condition face
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    function ielement_l(self,ind) result(ielem_bc)
        class(bc_patch_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: ind

        integer(ik) :: ielem_bc

        ielem_bc = self%ielement_l_%at(ind)

    end function ielement_l
    !******************************************************************************************







    !>  Return the face index of a boundary condition face
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    function iface(self,ind) result(iface_bc)
        class(bc_patch_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: ind

        integer(ik) :: iface_bc

        iface_bc = self%iface_%at(ind)

    end function iface
    !******************************************************************************************











    !>  For a face in the bc_patch, (iface_bc), add an element that it is coupled with.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_coupled_element(self,iface_bc,idomain_g, idomain_l, ielement_g, ielement_l, proc)
        class(bc_patch_t),  intent(inout)   :: self
        integer(ik),        intent(in)      :: iface_bc
        integer(ik),        intent(in)      :: idomain_g
        integer(ik),        intent(in)      :: idomain_l
        integer(ik),        intent(in)      :: ielement_g
        integer(ik),        intent(in)      :: ielement_l
        integer(ik),        intent(in)      :: proc


        call self%idomain_g_coupled(iface_bc)%push_back(idomain_g)
        call self%idomain_l_coupled(iface_bc)%push_back(idomain_l)
        call self%ielement_g_coupled(iface_bc)%push_back(ielement_g)
        call self%ielement_l_coupled(iface_bc)%push_back(ielement_l)
        call self%proc_coupled(iface_bc)%push_back(proc)

    end subroutine add_coupled_element
    !*******************************************************************************************






    !>  For a face in the bc_patch, (iface_bc), return the number of elements it is coupled with.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2017
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------
    function ncoupled_elements(self,iface_bc) result(res)
        class(bc_patch_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: iface_bc

        integer(ik) :: res

        res = self%ielement_g_coupled(iface_bc)%size()

    end function ncoupled_elements
    !********************************************************************************************












end module type_bc_patch
