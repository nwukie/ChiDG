module type_bc_patch
    use mod_kinds,      only: ik
    use type_ivector,   only: ivector_t
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    type, public :: bc_patch_t

        type(ivector_t)                 :: idomain_l
        type(ivector_t)                 :: ielement_l
        type(ivector_t)                 :: iface_l
        type(ivector_t), allocatable    :: coupled_elements(:)

    contains

        procedure   :: add_face
        procedure   :: nfaces

        ! Return indices for a given bc face
        procedure   :: idomain
        procedure   :: ielement
        procedure   :: iface

    end type bc_patch_t
    !******************************************************************************************






contains




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine add_face(self,idomain,ielement,iface)
        class(bc_patch_t),  intent(inout)   :: self
        integer(ik),        intent(in)      :: idomain
        integer(ik),        intent(in)      :: ielement
        integer(ik),        intent(in)      :: iface

        integer(ik)                     :: nfaces, ierr


        ! Add face to bc_patch list
        call self%idomain_l%push_back(idomain)
        call self%ielement_l%push_back(ielement)
        call self%iface_l%push_back(iface)


        ! Extend coupling storage. Need a vector for each bc_face
        nfaces = self%nfaces()
        if (allocated(self%coupled_elements)) deallocate(self%coupled_elements)
        allocate(self%coupled_elements(nfaces), stat=ierr)
        if (ierr /= 0) call AllocationError


    end subroutine add_face
    !**************************************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !--------------------------------------------------------------------------------------------------
    function nfaces(self) result(nfaces_bc)
        class(bc_patch_t),  intent(inout)   :: self

        integer(ik) :: nfaces_bc

        nfaces_bc = self%iface_l%size()

    end function nfaces
    !***************************************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    function idomain(self,ind) result(idom_bc)
        class(bc_patch_t),  intent(inout)   :: self
        integer(ik),        intent(in)      :: ind

        integer(ik) :: idom_bc

        idom_bc = self%idomain_l%at(ind)

    end function idomain
    !**************************************************************************************************






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/30/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    function ielement(self,ind) result(ielem_bc)
        class(bc_patch_t),  intent(inout)   :: self
        integer(ik),        intent(in)      :: ind

        integer(ik) :: ielem_bc

        ielem_bc = self%ielement_l%at(ind)

    end function ielement
    !**************************************************************************************************







    !>
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

        iface_bc = self%iface_l%at(ind)

    end function iface
    !**************************************************************************************************





end module type_bc_patch
