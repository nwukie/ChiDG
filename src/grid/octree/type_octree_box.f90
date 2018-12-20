!
! Octree box derived type
!
! @author Eric Wolf
! @date 08/25/2017
!

module type_octree_box
    use mod_kinds,                  only: ik, rk
    use mod_constants,              only: ZERO, ONE, TWO
    implicit none

    type, public :: octree_box_t


        ! If space_dim = 1, this is a binary tree.
        ! If space_dim = 2, this is a quadtree.
        ! If space_dim = 3, this is an octree.
        integer(ik)                 :: space_dim

        ! refinement_dir(i) = 1 if the tree is refined along direction i.
        ! sum(refinement_dir) = space_dim
        integer(ik)                 :: refinement_dir(3)

        ! Tree vertex indices for self, parent and children boxes.
        integer(ik)                 :: box_ID           ! Index of the box in the tree%vertices(:) array

        integer(ik)                 :: start_index      ! Start index of points in the box
        integer(ik)                 :: end_index        ! End index of points in the box
        integer(ik)                 :: num_points       ! Number of points in the box


        integer(ik)                 :: parent_ID        ! Array of the box indices whose support contains the present box
        !type(ivector_t)             :: child_ID      ! Array fo the v-indices contained in the support of the present box
        integer(ik)                 :: child_ID(8) = -1

        logical                     :: is_leaf = .true.

        ! Geometric information describing the spatial extent of the box.
        real(rk)                    :: center(3) ! (xctr,yctr,zctr) coordinates of the center of the box
        real(rk)                    :: extent(3) ! (xlen, ylen, zlen) extent (half side lengths) of the box

    contains

        procedure :: init

    end type octree_box_t

contains

    subroutine init(self,   corner_point_in, side_lengths_in, start_index, end_index, num_points, is_leaf)
        class(octree_box_t), intent(inout)    :: self
        real(rk),           intent(in)      :: corner_point_in(3), side_lengths_in(3)
        integer(ik),        intent(in)      :: start_index, end_index, num_points
        logical,            intent(in)      :: is_leaf
        

        
        !self%refinement_dir = refinement_dir_in
        !self%box_ID = box_ID_in
        self%center = corner_point_in
        self%extent = side_lengths_in
        self%start_index = start_index
        self%end_index   = end_index
        self%num_points = num_points
        self%is_leaf = is_leaf

        
    end subroutine init

end module type_octree_box
