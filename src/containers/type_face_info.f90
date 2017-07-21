module type_face_info
    use mod_kinds,  only: ik


    !> Container for locating a given face in the mesh data, via indices.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------
    type, public :: face_info_t

        integer(ik)     :: idomain_g
        integer(ik)     :: idomain_l
        integer(ik)     :: ielement_g
        integer(ik)     :: ielement_l
        integer(ik)     :: iface        !< Element-local face index

    end type face_info_t
    !**********************************************************************************


    interface face_info
        module procedure face_info_constructor
    end interface


contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/21/2017
    !!
    !----------------------------------------------------------------------------------
    function face_info_constructor(idomain_g, idomain_l, ielement_g, ielement_l, iface) result(face_info)
        integer(ik),    intent(in)  :: idomain_g
        integer(ik),    intent(in)  :: idomain_l
        integer(ik),    intent(in)  :: ielement_g
        integer(ik),    intent(in)  :: ielement_l
        integer(ik),    intent(in)  :: iface

        type(face_info_t)   :: face_info

        face_info%idomain_g  = idomain_g
        face_info%idomain_l  = idomain_l
        face_info%ielement_g = ielement_g
        face_info%ielement_l = ielement_l
        face_info%iface      = iface

    end function face_info_constructor
    !**********************************************************************************




end module type_face_info
