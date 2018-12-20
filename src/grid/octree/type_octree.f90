!
! Octree derived type supporting 3D radius neighbor search for RBF and other applications
!
! Based on: 
!    J. Behley, V. Steinhage, A.B. Cremers. Efficient Radius Neighbor Search in Three-dimensional Point Clouds,
!    Proc. of the IEEE International Conference on Robotics and Automation (ICRA), 2015
!
! @author Eric Wolf
! @date 08/25/2017
!
module type_octree
    use mod_kinds,                      only: ik, rk
    use mod_constants,                  only: ZERO, ONE, HALF, TWO
    use type_octree_box,                only: octree_box_t
    use type_box_vector,                only: box_vector_t
    use type_octree_parameters,         only: octree_parameters_t
    use type_ivector,                   only: ivector_t
    implicit none

    type, public :: octree_t

        type(octree_parameters_t)       :: params

        type(box_vector_t)              :: boxes

        integer(ik)                     :: root_box_ID


        real(rk),   allocatable         :: global_nodes(:, :)
        logical                         :: has_nodes = .false.
        
        integer(ik),    allocatable     :: successors(:)
    contains

        procedure                       :: init
        procedure                       :: build_octree_depth_first
        procedure                       :: create_octant
        procedure                       :: radius_search

    end type

contains

    subroutine init(self, bucket_size, min_extent, refine_dir, copy_points)
        class(octree_t),    intent(inout)                       :: self
        integer(ik),        intent(in),     optional            :: bucket_size
        real(rk),           intent(in),     optional            :: min_extent
        integer(ik),        intent(in),     optional            :: refine_dir(3)
        logical,            intent(in),     optional            :: copy_points

        call self%params%init(bucket_size, min_extent, refine_dir, copy_points)

    end subroutine init

    !> Depth-first octree construction
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    08/17/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine build_octree_depth_first(self, global_nodes)
        class(octree_t),    intent(inout)   :: self
        real(rk),           intent(in)      :: global_nodes(:,:)

        real(rk) :: center(3), extent(3)
        real(rk) :: coord_min(3), coord_max(3), bounding_box_len(3)

        integer(ik) :: ibox, ichild, current_level_ID, next_level_ID, current_box_ID, total_number_boxes, number_boxes, &
        current_level, child_box_ID, inode, level_counter, root_box_ID, start_index, end_index, num_points


        ! Copy the nodes if requested
        if (self%params%copy_points) then
            self%global_nodes = global_nodes
            self%has_nodes = .true.
        end if

        !----------------------------------------------------------------------------------------------------
        !
        !                   Create the root level (bounding box)
        !

        !
        ! Determine a bounding box for all nodes to construct the root box of the tree
        !
        coord_min(1) = minval(global_nodes(:,1))
        coord_min(2) = minval(global_nodes(:,2))
        coord_min(3) = minval(global_nodes(:,3))

        coord_max(1) = maxval(global_nodes(:,1))
        coord_max(2) = maxval(global_nodes(:,2))
        coord_max(3) = maxval(global_nodes(:,3))

        extent = HALF*(coord_max - coord_min)

        center = HALF*(coord_min + coord_max)
            
        num_points = size(global_nodes(:,1))
        start_index = 1
        end_index   = num_points

        allocate(self%successors(num_points))
        do inode = 1, num_points
            self%successors(inode) = inode+1
        end do



        !----------------------------------------------------------------------------------------------------
        !
        ! Now recursively subdivide boxes until the no box contains more than the max number of nodes per box
        !    


        call self%create_octant(center, extent, start_index, end_index, num_points, global_nodes, root_box_ID)
        self%root_box_ID = root_box_ID

    end subroutine build_octree_depth_first

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    08/17/2018 
    !!
    !--------------------------------------------------------------------------------
    recursive subroutine create_octant(self, center, extent, start_index, end_index, num_points, global_nodes, my_box_ID) 
        class(octree_t),    intent(inout)   :: self 
        real(rk),           intent(in)      :: center(3), extent(3)
        integer(ik),        intent(in)      :: num_points, start_index, end_index
        real(rk),           intent(in)      :: global_nodes(:,:)
        integer(ik),        intent(out)     :: my_box_ID

        
        integer(ik) :: box_ID, child_box_ID, idx, morton_code,last_child_box_ID, ichild, num_children, ipt
        integer(ik) :: child_starts(8), child_ends(8), child_sizes(8)
        type(octree_box_t)  :: temp_box
        logical     :: is_leaf, firsttime
        real(rk)    :: node(3), factor(2), child_center(3), child_extent(3)
        
        is_leaf = .true.
        !print *, 'a'
        call temp_box%init(center, extent, start_index, end_index, num_points, is_leaf)
        !print *, 'b'
        call self%boxes%push_back(temp_box, box_ID)
        my_box_ID = box_ID
        
        !print *, 'c'
        if (num_points > self%params%bucket_size) then
            ! Box is flagged for refinement
            self%boxes%data(box_ID)%is_leaf = .false.

            ! Loop over the points contained in current octant box and divide them into the child boxes
            idx = start_index
            !print *, 'd'
            child_sizes = 0
            child_starts = 1
            child_ends = 1
            do ipt = 1, num_points
                node = global_nodes(idx,:)

                ! Morton code trickery
                morton_code = 0
                if (node(1) > center(1)) morton_code = ior(morton_code, 1) 
                if (node(2) > center(2)) morton_code = ior(morton_code, 2) 
                if (node(3) > center(3)) morton_code = ior(morton_code, 4) 
                morton_code = morton_code + 1
                
                ! set child starts and update successors...
                if (child_sizes(morton_code) == 0) then
                    child_starts(morton_code) = idx
                else
                    self%successors(child_ends(morton_code)) = idx
                end if
                child_sizes(morton_code) = child_sizes(morton_code) + 1
                
                child_ends(morton_code) = idx
                idx = self%successors(idx)
            end do

            !print *, 'e'
            ! Loop over the non-empty child octants and create boxes
            child_extent = 0.5_rk*extent
            firsttime = .true.
            !last_child_idx = 1
            last_child_box_ID = 0
            factor = (/ -0.5_rk, 0.5_rk/)
            do ichild = 1, 8
                if (child_sizes(ichild) > 0) then

                    ! Get the child center using some Morton code trickery...
                    ! transfer(logical, int) converts the logical argument into integer
                    ! transfer(.false., 1) = 0, transfer(.true., 1) = 1
                    child_center(1) = center(1) + factor(transfer((iand(ichild-1,1)>0),1)+1)*extent(1)
                    child_center(2) = center(2) + factor(transfer((iand(ichild-1,2)>0),1)+1)*extent(2)
                    child_center(3) = center(3) + factor(transfer((iand(ichild-1,4)>0),1)+1)*extent(3)

                    ! Recursively call create_octant to create the child octant box
                    !print *, 'f'
                    call self%create_octant(child_center, child_extent, child_starts(ichild), child_ends(ichild), child_sizes(ichild), global_nodes, child_box_ID) 

                    !print *, 'g'
                    ! Store the child_box_ID in the current box
                    self%boxes%data(box_ID)%child_ID(ichild) = child_box_ID
                    !print *, 'g2'
                    if (firsttime) then
                    
                        ! Reset the start index of the current box to the start index of the first child
                        !print *, 'g3'
                        self%boxes%data(box_ID)%start_index = self%boxes%data(child_box_ID)%start_index
                        !print *, 'g4'
                    else
                      !successors_[octant->child[lastChildIdx]->end] =
                      !    octant->child[i]->start;  // we have to ensure that also the child ends link to the next child start.
                        !print *, 'g5'
                        self%successors(self%boxes%data(last_child_box_ID)%end_index) = self%boxes%data(child_box_ID)%start_index
                        !print *, 'g6'
                    
                    end if
                    !print *, 'h'
                    last_child_box_ID = child_box_ID
                    !lastChildIdx = i
                    !octant->end = octant->child[i]->end;
                    self%boxes%data(box_ID)%end_index = self%boxes%data(child_box_ID)%end_index
                    firsttime = .false.
                end if

            end do

        else
            
            !print *, 'leaf box'
        
        end if
        

    end subroutine create_octant


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    08/21/2018 
    !!
    !--------------------------------------------------------------------------------

    recursive subroutine radius_search(self, global_nodes, box_ID, query_point, query_radius, hit_list)
        class(octree_t),    intent(inout)   :: self
        real(rk),           intent(in)      :: global_nodes(:,:)
        integer(ik),        intent(in)      :: box_ID
        real(rk),           intent(in)      :: query_point(3), query_radius(3)
        type(ivector_t),    intent(inout)   :: hit_list

        real(rk) :: corner(3), center(3), extent(3), factor(2), point(3), child_center(3), child_extent(3)

        logical :: box_is_contained, corner_is_contained, point_included, child_overlaps
        integer(ik) :: icorner, idx, ipt, child_ID, ichild

        
        center = self%boxes%data(box_ID)%center
        extent = self%boxes%data(box_ID)%extent
        !print *, 'num_points'
        !print *, self%boxes%data(box_ID)%num_points

        ! Check if the current box is entirely contained in the searchoid
        ! If so, add all of the points in the box to the hit list
        factor = (/-1.0_rk, 1.0_rk/)
        box_is_contained = .true.
        do icorner = 1, 8
            corner(1) = center(1) + factor(transfer((iand(icorner-1,1)>0),1)+1)*extent(1)
            corner(2) = center(2) + factor(transfer((iand(icorner-1,2)>0),1)+1)*extent(2)
            corner(3) = center(3) + factor(transfer((iand(icorner-1,4)>0),1)+1)*extent(3)

            corner_is_contained = check_point_inclusion(corner, query_point, query_radius)
!            corner_is_contained = &
!            ((((corner(1)-query_point(1))/query_radius(1))**TWO+ &
!              ((corner(2)-query_point(2))/query_radius(2))**TWO+ &
!              ((corner(3)-query_point(3))/query_radius(3))**TWO)<=ONE)

            if (.not. (corner_is_contained)) box_is_contained = .false.
        end do
       
        if (box_is_contained) then
            ! Add all points to the hit list

            idx = self%boxes%data(box_ID)%start_index
            do ipt = 1, self%boxes%data(box_ID)%num_points
                call hit_list%push_back(idx)

                idx = self%successors(idx)

            end do


        
        ! Check if the current box is a leaf
        ! If so, individually check each point to see if it is contained in the searchoid
        ! and add contained points to the hit list
        ! Note: there should be just a small number points, at most the bucket size, in a leaf box

        else if (self%boxes%data(box_ID)%is_leaf) then

            ! Loop over all points in the present leaf box and test them for inclusion
            idx = self%boxes%data(box_ID)%start_index
            do ipt = 1, self%boxes%data(box_ID)%num_points
                point = global_nodes(idx,:)
                point_included = check_point_inclusion(point, query_point, query_radius) 
                if (point_included) call hit_list%push_back(idx)

                idx = self%successors(idx)

            end do

        ! Else, cycle over the child boxes of the current box,
        ! check for overlap with the searchoid, and if overlap is detected
        ! recursively launch this subroutine on the child box

        else
            do ichild = 1,8

                child_ID = self%boxes%data(box_ID)%child_ID(ichild)
                if (child_ID /= -1) then
                    ! Check overlap by putting bounding box around the searchoid
                    ! Then it's trivial to check whether the BB overlaps with the octree child box
                    child_center = self%boxes%data(child_ID)%center
                    child_extent = self%boxes%data(child_ID)%extent
                    child_overlaps = check_box_overlaps(child_center, child_extent, query_point, query_radius)

                    if (child_overlaps) call self%radius_search(global_nodes, child_ID, query_point, query_radius, hit_list)
                end if

            end do
        end if

    end subroutine radius_search

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    08/21/2018 
    !!
    !--------------------------------------------------------------------------------
    function check_point_inclusion(point, query_point, query_radius) result(point_included)
        real(rk), intent(in) :: point(3), query_point(3), query_radius(3)

        logical :: point_included

        point_included = &
            ((((point(1)-query_point(1))/query_radius(1))**TWO+ &
              ((point(2)-query_point(2))/query_radius(2))**TWO+ &
              ((point(3)-query_point(3))/query_radius(3))**TWO)<=ONE)


    end function check_point_inclusion


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    08/21/2018 
    !!
    !--------------------------------------------------------------------------------
    function check_box_overlaps(center, extent, query_point, query_radius) result(box_overlaps)
        real(rk)    :: center(3), extent(3), query_point(3), query_radius(3)
        logical     :: box_overlaps

        real(rk)    :: box_min(3), box_max(3), bbox_min(3), bbox_max(3)
        integer(ik) :: idir
        logical     :: bo(3)

        box_min = center - extent
        box_max = center + extent

        bbox_min = query_point - query_radius
        bbox_max = query_point + query_radius
        
        box_overlaps = .false.
        bo = .false.
        do idir = 1,3
            if (((bbox_min(idir)<=box_max(idir)) .and. (box_min(idir)<=bbox_min(idir))) .or. &
                ((box_min(idir)<=bbox_max(idir)) .and. (bbox_min(idir)<=box_min(idir)))) bo(idir) = .true.

        end do

        if ((bo(1)) .and. (bo(2)) .and. (bo(3))) box_overlaps = .true.

    end function check_box_overlaps
end module type_octree
