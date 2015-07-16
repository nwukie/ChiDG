module mod_grid_tools
    use mod_kinds,              only: rk,ik
    use mod_quadrature,         only: GQ
    use type_point,             only: point_t
    use type_expansion,         only: expansion_t
    use mod_element_mapping,    only: elem_map

    implicit none



contains


    !>  Compute discrete coordinate values, given the modes of a coordinate expansion
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  cmodes  Modal values for the coordinate expansion
    !!  @param[in]  igq     Integer for selecting the appropriate quadrature instance
    !!  @param[in]  pts     Array of computed coordinate points
    !-------------------------------------------------------------------------
    subroutine compute_discrete_coordinates(cmodes,igq,pts)
        type(expansion_t),              intent(in)      :: cmodes
        integer(ik),                    intent(in)      :: igq
        type(point_t),  dimension(:),   intent(inout)   :: pts


        pts(:)%c1_ = matmul(GQ(igq)%vol%val,cmodes%mat(:,1))
        pts(:)%c2_ = matmul(GQ(igq)%vol%val,cmodes%mat(:,2))
        pts(:)%c3_ = matmul(GQ(igq)%vol%val,cmodes%mat(:,3))

    end subroutine





    !>  Compute the modal representation of the point coordinates given an array of points
    !!  and an element mapping
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      pts     Array of points
    !!  @param[in]      imap    Integer for selecting the appropriate element mapping from 'elem_map'
    !!  @param[inout]   cmodes  Modal values for the coordinate expansion
    !-------------------------------------------------------------------------
    subroutine compute_modal_coordinates(pts,imap,cmodes)
        type(point_t),  dimension(:),   intent(in)    :: pts
        integer(ik),                    intent(in)    :: imap
        type(expansion_t),              intent(inout) :: cmodes

        if (size(elem_map(imap)%mat,1) /= size(pts)) stop "Error: compute_modal_coordinates -- mapping and point sizes do not match"

        cmodes%mat(:,1) = matmul(elem_map(imap)%mat,pts(:)%c1_)
        cmodes%mat(:,2) = matmul(elem_map(imap)%mat,pts(:)%c2_)
        cmodes%mat(:,3) = matmul(elem_map(imap)%mat,pts(:)%c3_)

    end subroutine
















end module mod_grid_tools
