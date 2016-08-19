module type_BR2_face_data
#include <messenger.h>
    use mod_kinds,          only: ik
    use mod_constants,      only: BR2_INTERIOR, BR2_EXTERIOR, INTERIOR, CHIMERA, BOUNDARY
    use type_mesh,          only: mesh_t
    use type_element_info,  only: element_info_t
    use type_chidgVector,   only: chidgVector_t
    use type_BR2_lift_diff, only: BR2_lift_diff_t
    implicit none





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !---------------------------------------------------------------------------------------
    type, public :: BR2_face_data_t

        type(BR2_lift_diff_t),   allocatable :: diff(:) ! lifting operators differentiated wrt 
                                                        ! dependent elements

    contains

        procedure   :: init
        procedure   :: update

    end type BR2_face_data_t
    !***************************************************************************************










contains








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self,mesh,elem_info,iface)
        class(BR2_face_data_t),     intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh(:)
        type(element_info_t),       intent(in)      :: elem_info
        integer(ik),                intent(in)      :: iface

        integer(ik) :: ndepend, ierr, ChiID, idiff
        logical     :: conforming_face, chimera_face, boundary_face, allocate_diff, reallocate_diff

        
        !
        ! Get number of dependent elements
        !
        conforming_face = (mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%ftype == INTERIOR)
        chimera_face    = (mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%ftype == CHIMERA )
        boundary_face   = (mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%ftype == BOUNDARY)


        if (conforming_face) then
            ! 1 Interior + 1 Exterior
            ndepend = 2

        else if (chimera_face) then
            ! 1 Interior + N chimera donors
            ChiID = mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%ChiID
            ndepend = 1 + mesh(elem_info%idomain_l)%chimera%recv%data(ChiID)%ndonors()

        else if (boundary_face) then
            ! 1 Interior
            ndepend = 1
        else
            call chidg_signal(FATAL,"BR2_lift%init: Invalid face type for computing BR2 lifting operators.")
        end if

        

        !
        ! Allocate storage for number of elements to be differentiated wrt
        !
        allocate_diff = (.not. allocated(self%diff))

        if (allocate_diff) then
            allocate(self%diff(ndepend), stat=ierr)
            if (ierr /= 0) call AllocationError

        else

            reallocate_diff = (ndepend /= size(self%diff))

            if (reallocate_diff) then
                deallocate(self%diff)
                allocate(self%diff(ndepend), stat=ierr)
                if (ierr /= 0) call AllocationError
            end if

        end if




        !
        ! Initialize storage for each differentiated element
        !
        do idiff = 1,size(self%diff)
            call self%diff(idiff)%init(mesh,elem_info,iface,idiff)
        end do


    end subroutine init
    !******************************************************************************************
















    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------
    subroutine update(self,mesh,elem_info,iface,q,BR2_TYPE)
        class(BR2_face_data_t),     intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh(:)
        type(element_info_t),       intent(in)      :: elem_info
        integer(ik),                intent(in)      :: iface
        type(chidgVector_t),        intent(in)      :: q
        integer(ik),                intent(in)      :: BR2_TYPE

        integer(ik) :: idiff

        !
        ! Compute the lifting operators differentiated with respect to each dependent element
        !
        do idiff = 1,size(self%diff)
            call self%diff(idiff)%update(mesh,elem_info,iface,idiff,q,BR2_TYPE)
        end do !idiff

    end subroutine update
    !**********************************************************************************************









end module type_BR2_face_data
