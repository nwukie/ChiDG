module type_BR2_face
    use mod_kinds,          only: ik
    use mod_constants,      only: BR2_INTERIOR, BR2_EXTERIOR, INTERIOR, CHIMERA
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

        logical :: conforming_face, chimera_face

    
        ! 
        ! Initialize/Reinitialize storage allocation and update interior lifting operators
        !
        call self%owner(1)%init(mesh,elem_info,iface)
        call self%owner(1)%update(mesh,elem_info,iface,q,BR2_INTERIOR)




        !
        ! Initialize/Reinitialize storage allocation and update exterior lifting operator
        ! if there is an exterior element
        !
        conforming_face = (mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%ftype == INTERIOR)
        chimera_face    = (mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%ftype == CHIMERA )
        if (conforming_face .or. chimera_face) then
            call self%owner(2)%init(mesh,elem_info,iface)
            call self%owner(2)%update(mesh,elem_info,iface,q,BR2_EXTERIOR)
        end if



    end subroutine update
    !****************************************************************************************





end module type_BR2_face
