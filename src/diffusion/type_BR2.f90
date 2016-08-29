module type_BR2
#include <messenger.h>
    use mod_kinds,          only: ik
    use mod_constants,      only: NFACES, BOUNDARY, DIAG, ME, NEIGHBOR, ZERO
    use type_mesh,          only: mesh_t
    use type_element_info,  only: element_info_t
    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t
    use type_chidgVector,   only: chidgVector_t
    use type_BR2_face,      only: BR2_face_t
    use mod_interpolate,    only: interpolate_face_autodiff, interpolate_element_autodiff
    use DNAD_D
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
        procedure   :: interpolate_lift_face
        procedure   :: interpolate_lift_element

        !procedure   :: interpolate_derivative_face
        !procedure   :: interpolate_derivative_element
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







    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    function interpolate_lift_element(self,mesh,elem_info,fcn_info,ieqn,interpolation_type) result(lift_gq)
        class(BR2_t),           intent(in)      :: self
        type(mesh_t),           intent(in)      :: mesh(:)
        type(element_info_t),   intent(in)      :: elem_info
        type(function_info_t),  intent(in)      :: fcn_info
        integer(ik),            intent(in)      :: ieqn
        character(len=*),       intent(in)      :: interpolation_type


        logical     :: valid(NFACES)
        integer(ik) :: lift_index, idiff_BR2, iface


        type(AD_D), dimension(mesh(elem_info%idomain_l)%elems(elem_info%ielement_l)%gq%vol%nnodes) :: &
                lift_gq


        select case (interpolation_type)
            case('ddx')
                lift_index = 1
            case('ddy')
                lift_index = 2
            case('ddz')
                lift_index = 3
            case default
                call chidg_signal(FATAL,"BR2%interpolate_lift_element: Invalid interpolation type. Options are 'ddx','ddy','ddz'.")
        end select


        if (fcn_info%idiff == DIAG) then
            idiff_BR2 = 1
        else
            idiff_BR2 = 1 + fcn_info%idepend
        end if

        associate ( val => mesh(elem_info%idomain_l)%elems(elem_info%ielement_l)%gq%vol%val )
        !
        ! Allocate lift_gq
        !
        do iface = 1,NFACES
            valid(iface) = (mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%ftype /= BOUNDARY)
            if (valid(iface)) then
                lift_gq = matmul(val, self%face(iface)%owner(1)%diff(idiff_BR2)%eqn(ieqn)%lift(:,lift_index))
                lift_gq = ZERO
                exit
            end if
        end do

        !
        ! Element Lift: R = r1 + r2 + r3 + r4 + r5 + r6
        !
        do iface = 1,NFACES
            valid(iface) = (mesh(elem_info%idomain_l)%faces(elem_info%ielement_l,iface)%ftype /= BOUNDARY)
            if (valid(iface)) then
                lift_gq = lift_gq + matmul(val, self%face(iface)%owner(1)%diff(idiff_BR2)%eqn(ieqn)%lift(:,lift_index))
            end if
        end do
        end associate



    end function interpolate_lift_element
    !************************************************************************************************









    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !------------------------------------------------------------------------------------------------------------------------
    function interpolate_lift_face(self,mesh,face_info,fcn_info,ieqn,interpolation_type,interpolation_source) result(lift_gq)
        class(BR2_t),           intent(in)      :: self
        type(mesh_t),           intent(in)      :: mesh(:)
        type(face_info_t),      intent(in)      :: face_info
        type(function_info_t),  intent(in)      :: fcn_info
        integer(ik),            intent(in)      :: ieqn
        character(len=*),       intent(in)      :: interpolation_type
        integer(ik),            intent(in)      :: interpolation_source


        logical     :: valid(NFACES)
        integer(ik) :: lift_index, idiff_BR2, iface_n


        type(AD_D), dimension(mesh(face_info%idomain_l)%elems(face_info%ielement_l)%gq%face%nnodes) :: &
                lift_gq


        select case (interpolation_type)
            case('ddx')
                lift_index = 1
            case('ddy')
                lift_index = 2
            case('ddz')
                lift_index = 3
            case default
                call chidg_signal(FATAL,"BR2%interpolate_lift_element: Invalid interpolation type. Options are 'ddx','ddy','ddz'.")
        end select


        if (fcn_info%idiff == DIAG) then
            idiff_BR2 = 1
        else
            idiff_BR2 = 1 + fcn_info%idepend
        end if


        !
        ! Face Lift: r_i
        !
        if (interpolation_source == ME) then
            associate ( val => mesh(face_info%idomain_l)%elems(face_info%ielement_l)%gq%face%val(:,:,face_info%iface) )
                lift_gq = matmul(val, self%face(face_info%iface)%owner(1)%diff(idiff_BR2)%eqn(ieqn)%lift(:,lift_index) )
            end associate



        else if (interpolation_source == NEIGHBOR) then
            iface_n = mesh(face_info%idomain_l)%faces(face_info%ielement_l,face_info%iface)%get_neighbor_face()

            associate ( val => mesh(face_info%idomain_l)%elems(face_info%ielement_l)%gq%face%val(:,:,iface_n) )
                lift_gq = matmul(val, self%face(face_info%iface)%owner(2)%diff(idiff_BR2)%eqn(ieqn)%lift(:,lift_index) )
            end associate


        end if





    end function interpolate_lift_face
    !**********************************************************************************************************









!    !>
!    !!
!    !!
!    !!
!    !!
!    !!
!    !!
!    !---------------------------------------------------------------------------------------------------------
!    function interpolate_derivative_face()
!
!
!
!        dg_deriv = interpolate_face_autodiff()
!
!        lift = interpolate_lift_face
!
!        
!        br2_deriv = dg_deriv + lift
!
!
!
!    end function interpolate_derivative_face
!    !*********************************************************************************************************





end module type_BR2
