module type_BR2
    use mod_kinds,          only: ik
    use mod_constants,      only: NFACES
    use type_mesh,          only: mesh_t
    use type_element_info,  only: element_info_t
    use type_chidgVector,   only: chidgVector_t
    use type_BR2_face,      only: BR2_face_t
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !-------------------------------------------------------------------------------
    type, public :: BR2_t

        type(BR2_face_t)    :: face(NFACES)

    contains

        procedure   :: update

    end type BR2_t
    !*******************************************************************************






contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine update(self,mesh,elem_info,q)
        class(BR2_t),           intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh(:)
        type(element_info_t),   intent(in)      :: elem_info
        type(chidgVector_t),    intent(in)      :: q

        integer(ik) :: iface

        ! Compute lifting operators along each face
        do iface = 1,size(self%face)
            call self%face(iface)%update(mesh,elem_info,iface,q)
        end do !iface

    end subroutine update
    !***********************************************************************************************







end module type_BR2
