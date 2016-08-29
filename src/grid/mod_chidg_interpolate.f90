module mod_chidg_interpolate
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use DNAD_D

    use type_mesh,          only: mesh_t
    use type_chidgVector,   only: chidgVector_t
    use type_element_info,  only: element_info_t
    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t
    use mod_interpolate,    only: interpolate_element_autodiff, interpolate_face_autodiff
    implicit none


    interface interpolate
        module procedure    chidg_interpolate_element, chidg_interpolate_face
    end interface


contains




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------------------
    function chidg_interpolate_element(mesh,q,elem_info,fcn_info,ieqn,interpolation_type) result(var_gq)
        type(mesh_t),           intent(in)      :: mesh(:)
        type(chidgVector_t),    intent(in)      :: q
        type(element_info_t),   intent(in)      :: elem_info
        type(function_info_t),  intent(in)      :: fcn_info
        integer(ik),            intent(in)      :: ieqn
        character(len=*),       intent(in)      :: interpolation_type

        type(AD_D), dimension(mesh(elem_info%idomain_l)%elems(elem_info%ielement_l)%gq%vol%nnodes)  :: &
            var_gq, diff, lift

        if (interpolation_type == 'value') then

            var_gq = interpolate_element_autodiff(mesh,q,elem_info,fcn_info,ieqn,interpolation_type)

        else if (interpolation_type == 'ddx' .or. &
                 interpolation_type == 'ddy' .or. & 
                 interpolation_type == 'ddz') then

            diff = interpolate_element_autodiff(mesh,q,elem_info,fcn_info,ieqn,interpolation_type)
!            lift = BR2%interpolate_lift_element(mesh,elem_info,fcn_info,ieqn,interpolation_type)

            var_gq = diff + lift

        else
            call chidg_signal(FATAL,"interpolate: Invalid interpolation type. Options are 'value', 'ddx', 'ddy', 'ddz'")
        end if



    end function chidg_interpolate_element
    !********************************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/19/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------------------
    function chidg_interpolate_face(mesh,q,face_info,fcn_info,ieqn,interpolation_type,interpolation_source) result(var_gq)
        type(mesh_t),           intent(in)      :: mesh(:)
        type(chidgVector_t),    intent(in)      :: q
        type(face_info_t),      intent(in)      :: face_info
        type(function_info_t),  intent(in)      :: fcn_info
        integer(ik),            intent(in)      :: ieqn
        character(len=*),       intent(in)      :: interpolation_type
        integer(ik),            intent(in)      :: interpolation_source

        type(AD_D), dimension(mesh(face_info%idomain_l)%elems(face_info%ielement_l)%gq%face%nnodes) :: &
                var_gq, diff, lift

        if (interpolation_type == 'value') then

            var_gq = interpolate_face_autodiff(mesh,q,face_info,fcn_info,ieqn,interpolation_type,interpolation_source)

        else if (interpolation_type == 'ddx' .or. &
                 interpolation_type == 'ddy' .or. & 
                 interpolation_type == 'ddz') then

            diff = interpolate_face_autodiff(mesh,q,face_info,fcn_info,ieqn,interpolation_type,interpolation_source)
!            lift = BR2%interpolate_lift_face(mesh,face_info,fcn_info,ieqn,interpolation_type,interpolation_source)

            var_gq = diff + lift

            call chidg_signal(FATAL,"interpolate: Invalid interpolation type. Options are 'value', 'ddx', 'ddy', 'ddz'")
        end if



    end function chidg_interpolate_face
    !********************************************************************************************************











end module mod_chidg_interpolate
