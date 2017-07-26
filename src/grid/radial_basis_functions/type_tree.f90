module type_tree
    use mod_kinds, only: ik, rk
    implicit none

    type, public :: tree_t

        type(tree_root_t) :: root
        type(tree_vertex_t), allocatable :: vertices(:)

    contains

        procedure :: init_tree

        procedure   :: construct_tree
        procedure   :: compute_coefficients
        procedure   :: evalutate_rbf

    end type tree_t


contains

    subroutine construct_tree(self, nodes, base_nodes_indices, support_radius, rbf_name)
        type(tree_t), intent(inout) :: self
        real(rk),       intent(inout)   :: nodes(:)
        integer(ik),    intent(in)      :: base_nodes_indices
        real(rk), intent(in)                            :: support_radius
        character(:), intent(in)                        :: rbf_name


        integer(ik) :: inode, inode_base
        ! Process base nodes

        do inode_base = 1, size(base_nodes_indices)

            inode = base_node_indices(inode_base)

            self%root%add_base_node(nodes(inode),support_radius,rbf_name)
        end do

        ! Add remaining surface points in a refinement
        num_active_points = size(base_nodes_indices)
        active_nodes_indices(1:num_active_nodes) = base_nodes_indices

        do while (num_active_points<size(nodes))
            !Compute separations


            !Find maximum separation

            !Add to the active point list and consuct tree
            
            !If the new refinement point is contained only in the base point
            ! and not in any other refinement point support, then it is a Level 1 vertex.

            !If the new refinement point is contained in the support of a Level k vertex,
            ! then it is a Level k+1 vertex.

            !How do we know which case occurs?
            ! Look at the nearest points.
            ! If the nearest points are all base points, then Level 1.
            ! Else if a nearest point is a Level k vertex, then it is a Level k+1 vertex.


        end do

    end subroutine


end module type_tree
