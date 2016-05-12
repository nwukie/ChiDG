module mod_grid_tools
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_quadrature,         only: GQ
    use type_point,             only: point_t
    use type_densevector,       only: densevector_t
    use mod_grid,               only: get_element_mapping
    implicit none





contains




    !>  Compute discrete coordinate values, given the modes of a coordinate expansion
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  cmodes  Modal values for the coordinate expansion
    !!  @param[in]  igq     Integer for selecting the appropriate quadrature instance
    !!  @param[in]  pts     Array of computed coordinate points
    !!
    !---------------------------------------------------------------------------------------------------------
    subroutine compute_discrete_coordinates(cmodes,igq,pts)
        type(densevector_t),            intent(in)      :: cmodes
        integer(ik),                    intent(in)      :: igq
        type(point_t),  dimension(:),   intent(inout)   :: pts


        pts(:)%c1_ = matmul(GQ(igq)%vol%val,cmodes%getvar(1))
        pts(:)%c2_ = matmul(GQ(igq)%vol%val,cmodes%getvar(2))
        pts(:)%c3_ = matmul(GQ(igq)%vol%val,cmodes%getvar(3))

    end subroutine compute_discrete_coordinates
    !*********************************************************************************************************









    !>  Compute the modal representation of the point coordinates given an array of points
    !!  and an element mapping
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]      pts     Array of points
    !!  @param[in]      imap    Integer for selecting the appropriate element mapping from 'elem_map'
    !!  @param[inout]   cmodes  Modal values for the coordinate expansion
    !---------------------------------------------------------------------------------------------------------
    subroutine compute_modal_coordinates(spacedim,pts,imap,cmodes)
        integer(ik),                    intent(in)      :: spacedim
        type(point_t),  dimension(:),   intent(in)      :: pts
        integer(ik),                    intent(in)      :: imap
        type(densevector_t),            intent(inout)   :: cmodes

        real(rk), dimension(size(pts))  :: xmodes, ymodes, zmodes
        real(rk),   allocatable         :: element_mapping(:,:)



        element_mapping = get_element_mapping(spacedim,imap)
        if (size(element_mapping,1) /= size(pts) ) call chidg_signal(FATAL, "compute_modal_coordinate: mapping and point sizes do not match.")


        xmodes = matmul(element_mapping,pts(:)%c1_)
        ymodes = matmul(element_mapping,pts(:)%c2_)
        zmodes = matmul(element_mapping,pts(:)%c3_)



        call cmodes%setvar(1,xmodes)
        call cmodes%setvar(2,ymodes)
        call cmodes%setvar(3,zmodes)

    end subroutine compute_modal_coordinates
    !**********************************************************************************************************








end module mod_grid_tools
