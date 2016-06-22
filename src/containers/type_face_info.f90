module type_face_info
    use mod_kinds,  only: ik
    use type_seed,  only: seed_t



    !> Container for locating a given face in the mesh data, via indices.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------
    type, public :: face_info_t

!        integer(ik)     :: idomain      !< Domain index
!        integer(ik)     :: ielement     !< Domain-local element index
        integer(ik)     :: idomain_g
        integer(ik)     :: idomain_l
        integer(ik)     :: ielement_g
        integer(ik)     :: ielement_l
        integer(ik)     :: iface        !< Element-local face index
        type(seed_t)    :: seed         !< Indices of the element being linearized

    end type face_info_t
    !**********************************************************************************





end module type_face_info
