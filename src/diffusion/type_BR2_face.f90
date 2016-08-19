module type_BR2_face
    use mod_kinds,          only: ik
    use mod_constants,      only: BR2_INTERIOR, BR2_EXTERIOR
    use type_mesh,          only: mesh_t
    use type_element_info,  only: element_info_t
    use type_chidgVector,   only: chidgVector_t
    use type_BR2_face_data, only: BR2_face_data_t
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !----------------------------------------------------------------------------------------
    type, public :: BR2_face_t

        type(BR2_face_data_t)   :: owner(2)     ! 1 = BR2_INTERIOR, 2 = BR2_EXTERIOR

    contains

        procedure   :: update

    end type BR2_face_t
    !****************************************************************************************







contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine update(self,mesh,elem_info,iface,q)
        class(BR2_face_t),      intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh(:)
        type(element_info_t),   intent(in)      :: elem_info
        integer(ik),            intent(in)      :: iface
        type(chidgVector_t),    intent(in)      :: q


    
        ! 
        ! Initialize/Reinitialize storage allocation if necessary
        !
        call self%owner(1)%init(mesh,elem_info,iface)
        call self%owner(2)%init(mesh,elem_info,iface)



        !
        ! Update face data
        !
        call self%owner(1)%update(mesh,elem_info,iface,q,BR2_INTERIOR)
        call self%owner(2)%update(mesh,elem_info,iface,q,BR2_EXTERIOR)



    end subroutine update
    !****************************************************************************************





end module type_BR2_face
