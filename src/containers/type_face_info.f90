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

    contains

        procedure   :: init

    end type face_info_t
    !**********************************************************************************



contains




    !>  Initialize face info data.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/20/2016
    !!
    !----------------------------------------------------------------------------------
    subroutine init(self,idomain_g,idomain_l,ielement_g,ielement_l,iface)
        class(face_info_t), intent(inout)   :: self
        integer(ik),        intent(in)      :: idomain_g
        integer(ik),        intent(in)      :: idomain_l
        integer(ik),        intent(in)      :: ielement_g
        integer(ik),        intent(in)      :: ielement_l
        integer(ik),        intent(in)      :: iface


        self%idomain_g  = idomain_g
        self%idomain_l  = idomain_l
        self%ielement_g = ielement_g
        self%ielement_l = ielement_l
        self%iface      = iface

    end subroutine init
    !**********************************************************************************


end module type_face_info
