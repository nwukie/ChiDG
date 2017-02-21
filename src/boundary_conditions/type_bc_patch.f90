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
    !!
    !!
    !-----------------------------------------------------------------------------------------
    type, public :: bc_patch_t

        type(ivector_t)                 :: idomain_l_
        type(ivector_t)                 :: ielement_l_
        type(ivector_t)                 :: iface_
        type(ivector_t), allocatable    :: coupled_elements(:)

    contains

        procedure   :: add_face
        procedure   :: nfaces

        ! Return indices for a given bc face
        procedure   :: idomain_l
        procedure   :: ielement_l
        procedure   :: iface

    end type bc_patch_t
    !******************************************************************************************






contains




    !>  Add a face to the boundary condition patch.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    function add_face(self,idomain,ielement,iface) result(iface_bc)
        class(bc_patch_t),  intent(inout)   :: self
        integer(ik),        intent(in)      :: idomain
        integer(ik),        intent(in)      :: ielement
        integer(ik),        intent(in)      :: iface

        integer(ik)                     :: nfaces, ierr, iface_bc

        !
        ! Add face to bc_patch list
        !
        call self%idomain_l_%push_back(idomain)
        call self%ielement_l_%push_back(ielement)
        call self%iface_%push_back(iface)


        !
        ! Get location of face in bc_patch
        !
        iface_bc = self%iface_%size()


        !
        ! Extend coupling storage. Need a vector for each bc_face
        !
        nfaces = self%nfaces()
        if (allocated(self%coupled_elements)) deallocate(self%coupled_elements)
        allocate(self%coupled_elements(nfaces), stat=ierr)
        if (ierr /= 0) call AllocationError


    end function add_face
    !**************************************************************************************************






    !>  Return the number of faces in the boundary condition patch
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !--------------------------------------------------------------------------------------------------
    function nfaces(self) result(nfaces_bc)
        class(bc_patch_t),  intent(inout)   :: self

        integer(ik) :: nfaces_bc

        nfaces_bc = self%iface_%size()

    end function nfaces
    !***************************************************************************************************








    !>  Return the local domain index of a boundary condition face
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    function idomain_l(self,ind) result(idom_bc)
        class(bc_patch_t),  intent(inout)   :: self
        integer(ik),        intent(in)      :: ind

        integer(ik) :: idom_bc

        idom_bc = self%idomain_l_%at(ind)

    end function idomain_l
    !**************************************************************************************************






    !>  Return the local element index of a boundary condition face
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    function ielement_l(self,ind) result(ielem_bc)
        class(bc_patch_t),  intent(inout)   :: self
        integer(ik),        intent(in)      :: ind

        integer(ik) :: ielem_bc

        ielem_bc = self%ielement_l_%at(ind)

    end function ielement_l
    !**************************************************************************************************







    !>  Return the face index of a boundary condition face
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    function iface(self,ind) result(iface_bc)
        class(bc_patch_t),  intent(inout)   :: self
        integer(ik),        intent(in)      :: ind

        integer(ik) :: iface_bc

        iface_bc = self%iface_%at(ind)

    end function iface
    !**************************************************************************************************





end module type_bc_patch
