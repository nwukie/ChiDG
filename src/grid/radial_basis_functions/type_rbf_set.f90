module type_rbf_set
    implicit none

    !> The rbf_set_t derived type provides a general container for a set of 
    !! compactly-supported RBFs and tools to allow for their rapid evaluation
    !! using an octree-based radius search algorithm.
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    08/31/2018 
    !!
    !--------------------------------------------------------------------------------
    type, public :: rbf_set_t
        
        integer(ik) :: rbf_set_ID
        
        
        real(rk), allocatable :: center(:,:), radius(:,:) !(nrbf, 3)
        real(rk), allocatable :: coefficients(:,:) ! (nrbf, ndata)


    end type rbf_set_t


contains

    !> For RBF-pull operations (e.g., RBF-smoothing of AV), we do a pre-processing step 
    !! where the RBFs register themselves with the elements they touch. Later, when looping
    !! over elements, the elements know which RBFs to evaluate.
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    08/31/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine register_rbfs_with_nodes(self, octree)
        class(rbf_set_t), intent(inout) :: self
        type(octree_t),     intent(inout) :: octree

        type(ivector_t) :: hit_list
        type(element_info_vector_t) :: elem_info_vec

        do irbf = 1, size(self%center(:,1))

            call hit_list%clear()
            call octree%radius_search(octree%global_nodes, octree%root_box_ID, self%center(irbf, :), self%radius(irbf,:), hit_list)

            do inode = 1, size(hit_list)
                
                ! Get the elem_info for all elements associated with the node in the global_node list
                elem_info_vec = get_element_info(hit_list%at(inode))


                ! Loop over all the elements assocated with the node and register the current RBF 
                do ielem = 1, size(elem_info_vec)
                    call register_rbf_with_elem(elem_info%at(ielem), rbf_set_ID, irbf)
                end do

            end do

        end do

    end subroutine register_rbfs_with_nodes

end module type_rbf_set
