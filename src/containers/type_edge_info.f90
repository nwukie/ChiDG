module type_edge_info
    use mod_kinds,  only: ik


    !>  Container for locating a given edge in the mesh data, via indices.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/6/2017
    !!
    !----------------------------------------------------------------------------------
    type, public :: edge_info_t

        integer(ik)     :: idomain_g
        integer(ik)     :: idomain_l
        integer(ik)     :: ielement_g
        integer(ik)     :: ielement_l
        integer(ik)     :: iedge

    end type edge_info_t
    !**********************************************************************************




end module type_edge_info
