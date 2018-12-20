module type_rbf_interpolation
#include <messenger.h>
    use mod_kinds,                              only: rk, ik
    use mod_constants,                          only: ZERO, ONE, TWO
    use mod_inv,                                only: inv
    use mod_radial_basis_function,              only: create_radial_basis_function
    use mod_rbf_tools
    use mod_k_means_clustering
    use type_radial_basis_function,             only: radial_basis_function_t
    use type_rbf_info,                          only: rbf_info_t
    use type_rbf_node,                          only: rbf_node_t
    use type_rbf_node_patch,                    only: rbf_node_patch_t
    use type_rbf_node_vector,                   only: rbf_node_vector_t
    use type_rbf_interpolation_matrix,          only: rbf_interpolation_matrix_t
    use type_ivector,                           only: ivector_t
    use type_svector,                           only: svector_t
    use type_mesh,                              only: mesh_t
    implicit none
    
    type :: rbf_interpolation_t

        type(rbf_info_t)                                    :: rbf_info

        class(radial_basis_function_t), allocatable         :: rbf

        ! For the multiscale RBFI, we might want two different kinds of RBFs
        ! for the base nodes and explicit nodes. We might want to use
        ! global RBFs (e.g. thin-plate splines) for the base nodes, 
        ! while the use of compactly supported RBFs for the explicit nodes seems to be mandatory.
        class(radial_basis_function_t), allocatable         :: rbf_explicit
    
        ! RBF node patch, with patch node ordering.
        type(rbf_node_patch_t)                              :: rbf_node_patch

        ! RBF node lists. Ordered according to RBF system node ordering.
        type(rbf_node_vector_t)                             :: source_nodes
        type(rbf_node_vector_t)                             :: source_nodes_base
        type(rbf_node_vector_t)                             :: source_nodes_explicit

        ! RBF interpolation matrix
        type(rbf_interpolation_matrix_t)                    :: interpolation_matrix      

    contains

        !procedure           :: init_from_patch 
        procedure           :: init_from_array
        procedure           :: init_from_array_std
        procedure           :: add_rbf
        procedure           :: add_rbf_explicit
        !procedure           :: assemble_node_patch
        procedure           :: assemble_rbf_patch
        procedure           :: assemble_node_patch_from_array
        procedure           :: construct_rbf_interpolation
        procedure           :: construct_rbf_interpolation_std
        procedure           :: solve
        procedure           :: evaluate_rbf
        procedure           :: evaluate_rbf_grad
        procedure           :: evaluate
        procedure           :: evaluate_grad
        procedure           :: evaluate_ms
        procedure           :: evaluate_ms_grad

    end type rbf_interpolation_t

contains

!    !>  Provides an interface to initialize an RBF interpolation from a specified 
!    !!  BC patch_group.
!    !!  
!    !! 
!    !!  @author Eric Wolf
!    !!  @date 10/19/2017
!    !--------------------------------------------------------------------------------
!    subroutine init_from_patch(self, rbfstring, rbfstring_explicit, data, patch_group, nnodes_base, radius_base)
!        class(rbf_interpolation_t),             intent(inout)   :: self
!        character(*),                           intent(in)      :: rbfstring
!        character(*),                           intent(in)      :: rbfstring_explicit
!        type(chidg_data_t),                     intent(inout)   :: data
!        character(*),                           intent(in)      :: patch_group
!        integer(ik),                            intent(in)      :: nnodes_base
!        real(rk),                               intent(in)      :: radius_base(3)
!
!        call self%add_rbf(rbfstring)
!        call self%add_rbf_explicit(rbfstring_explicit)
!        call self%assemble_node_patch(data, patch_group)
!        call self%construct_rbf_interpolation(nnodes_base, radius_base)
!
!    end subroutine init_from_patch
!    !********************************************************************************
!
    !>  Provides an interface to initialize an RBF interpolation from a specified 
    !!  BC patch_group.
    !!  
    !! 
    !!  @author Eric Wolf
    !!  @date 10/19/2017
    !--------------------------------------------------------------------------------
    subroutine init_from_array(self, rbfstring, rbfstring_explicit, source_nodes_in, nnodes_base, radius_base)
        class(rbf_interpolation_t),             intent(inout)   :: self
        character(*),                           intent(in)      :: rbfstring
        character(*),                           intent(in)      :: rbfstring_explicit
        real(rk),                               intent(in)      :: source_nodes_in(:,:)
        integer(ik),                            intent(in)      :: nnodes_base
        real(rk),                               intent(in)      :: radius_base(3)

        call self%add_rbf(rbfstring)
        call self%add_rbf_explicit(rbfstring_explicit)
        call self%assemble_node_patch_from_array(source_nodes_in)
        call self%construct_rbf_interpolation(nnodes_base, radius_base)

    end subroutine init_from_array
    !********************************************************************************

    !>  Provides an interface to initialize an RBF interpolation from a specified 
    !!  BC patch_group.
    !!  
    !! 
    !!  @author Eric Wolf
    !!  @date 10/19/2017
    !--------------------------------------------------------------------------------
    subroutine init_from_array_std(self, rbfstring, source_nodes_in, radius_base)
        class(rbf_interpolation_t),             intent(inout)   :: self
        character(*),                           intent(in)      :: rbfstring
        real(rk),                               intent(in)      :: source_nodes_in(:,:)
        real(rk),                               intent(in)      :: radius_base(3)

        call self%add_rbf(rbfstring)
        call self%assemble_node_patch_from_array(source_nodes_in)
        call self%construct_rbf_interpolation_std(radius_base)

    end subroutine init_from_array_std
    !********************************************************************************



    !>  Adds a rbf instance to self, choosing the type based on input string rbfstring
    !!  
    !! 
    !!  @author Eric Wolf
    !!  @date 10/19/2017
    !--------------------------------------------------------------------------------
    subroutine add_rbf(self, rbfstring)
        class(rbf_interpolation_t),             intent(inout)   :: self
        character(*),                           intent(in)      :: rbfstring

        class(radial_basis_function_t), allocatable                :: rbf
            
        integer(ik)     :: ierr
        !call self%set_rbf_name(rbfstring)
        call create_radial_basis_function(rbf, rbfstring)
        if (allocated(self%rbf)) deallocate(self%rbf)
        allocate(self%rbf, source = rbf, stat=ierr)
        if (ierr /= 0) call AllocationError

    end subroutine
    !********************************************************************************

    !>  Adds a rbf instance to self, choosing the type based on input string rbfstring
    !!  
    !! 
    !!  @author Eric Wolf
    !!  @date 10/19/2017
    !--------------------------------------------------------------------------------
    subroutine add_rbf_explicit(self, rbfstring)
        class(rbf_interpolation_t),             intent(inout)   :: self
        character(*),                           intent(in)      :: rbfstring

        class(radial_basis_function_t), allocatable                :: rbf
            
        integer(ik)     :: ierr
        !call self%set_rbf_name(rbfstring)
        call create_radial_basis_function(rbf, rbfstring)
        if (allocated(self%rbf_explicit)) deallocate(self%rbf_explicit)
        allocate(self%rbf_explicit, source = rbf, stat=ierr)
        if (ierr /= 0) call AllocationError

    end subroutine add_rbf_explicit
    !********************************************************************************


    !>  Assembles the rbf_node_patch from a BC patch_group.
    !!  
    !! 
    !!  @author Eric Wolf
    !!  @date 10/19/2017
    !--------------------------------------------------------------------------------
   subroutine assemble_rbf_patch(self,mesh,patch_group)
        class(rbf_interpolation_t),         intent(inout)       :: self
        type(mesh_t),                 intent(inout)       :: mesh 
        type(svector_t),                       intent(inout)          :: patch_group

        call self%rbf_node_patch%assemble_rbf_patch(mesh,patch_group)

    end subroutine assemble_rbf_patch
    !********************************************************************************


!    !>  Assembles the rbf_node_patch from a BC patch_group.
!    !!  
!    !! 
!    !!  @author Eric Wolf
!    !!  @date 10/19/2017
!    !--------------------------------------------------------------------------------
!   subroutine assemble_node_patch(self,data,patch_group)
!        class(rbf_interpolation_t),         intent(inout)       :: self
!        type(chidg_data_t),                 intent(inout)       :: data
!        character(*),                       intent(in)          :: patch_group
!
!        call self%rbf_node_patch%assemble(data,patch_group)
!
!    end subroutine assemble_node_patch
!    !********************************************************************************

    !>  Assembles self%rbf_node_patch from a given array of nodes.
    !!  
    !! 
    !!  @author Eric Wolf
    !!  @date 10/19/2017
    !--------------------------------------------------------------------------------
    subroutine assemble_node_patch_from_array(self,source_nodes_in) 
        class(rbf_interpolation_t),         intent(inout)       :: self
        real(rk),                           intent(in)          :: source_nodes_in(:,:)

        call self%rbf_node_patch%assemble_from_array(source_nodes_in)

    end subroutine assemble_node_patch_from_array
    !********************************************************************************

    !>  Constructs the compactly supported RBF interpolation according to 
    !!  the multiscale formulation of
    !!      Kedward, Allen, Rendall (2017)
    !!  
    !!  The source nodes are divided into a small set of base nodes,
    !!  numbering nnodes_base, and a much larger set of explicit nodes. The names
    !!  'base' and 'explicit' pertain to the construction and solution of the 
    !!  linear system used to obtain RBF coefficients from nodal values.
    !!
    !!  The base nodes give a dense RBF submatrix A, whose inverse Ainv we obtain
    !!  with a direct solver and store for later repeated use, while the 
    !!  explicit nodes give a dense block B in connection to the base nodes
    !!  and a lower triangular RBF submatrix C in connection with explicit nodes. 
    !!  The full RBF interpolation matrix has the form,
    !!  
    !!  M = 
    !!      |A 0|
    !!      |B C|.
    !!
    !!  We solve the linear system M*x=f arising in RBF interpolation with 
    !!  base and explicit nodal data,
    !!
    !!  f =
    !!      |f_b|
    !!      |f_e|,
    !!
    !!  and RBF coefficient solution vector,
    !!
    !!  x =
    !!      |x_b|
    !!      |x_e|,
    !!
    !!  by applying the stored Ainv,
    !!
    !!  x_b = Ainv*f_b
    !!
    !!  and then applying forward substitution to solve the lower-triangular system
    !!
    !!  C*x_e = f_e-B*x_b.
    !!
    !!  This linear system solution method costs O(N_b^3 + N_b*N_e + N_e) operations,
    !!  and we will generally take N_b << N_e to be a small fixed number.
    !!
    !!  We use k-means clustering with k=N_b to automatically select our base node set.
    !!
    !! 
    !!  @author Eric Wolf
    !!  @date 10/19/2017
    !--------------------------------------------------------------------------------
    subroutine construct_rbf_interpolation(self, nnodes_base, radius_base)
        class(rbf_interpolation_t),         intent(inout)   :: self
        integer(ik),                        intent(in), optional :: nnodes_base
        real(rk),                           intent(in), optional :: radius_base(3)

        type(ivector_t) :: active_nodes, inactive_nodes

        real(rk) :: radius_new(3), target_node(3), source_node(3), radius_rbf(3)
        real(rk) :: dist_max, min_dist, current_dist
        integer(ik) :: nnodes_total, nnodes_explicit
        integer(ik) :: inode, inode_active, inode_transfer, inode_ag, inode_ig, num_active, num_inactive, inode_source, inode_target

        type(rbf_node_t) :: rbf_node_temp

        integer(ik), allocatable :: base_node_indices(:)

        ! We may want to implement a NON-multiscale version to have a basis for comparision.
        if ((present(nnodes_base)) .and. (present(radius_base))) then

            ! Select the nnodes_base number of base nodes using k-means clustering
            nnodes_total = size(self%rbf_node_patch%nodes, 1)
            nnodes_explicit = nnodes_total-nnodes_base
            base_node_indices = get_k_means_node_indices(self%rbf_node_patch%nodes, nnodes_base)

            
            ! Store the base node informationinto the relevant data strucutres
            ! These comprise the initial active node set
            do inode = 1, nnodes_base

                call rbf_node_temp%set_node(self%rbf_node_patch%nodes(base_node_indices(inode),:), radius_base)
                call active_nodes%push_back(base_node_indices(inode))
                call self%source_nodes_base%push_back(rbf_node_temp)
                call self%source_nodes%push_back(rbf_node_temp)
                self%rbf_node_patch%patch_to_rbf_index(base_node_indices(inode)) = inode
                self%rbf_node_patch%rbf_to_patch_index(inode) = base_node_indices(inode)

            end do

            ! The remaining nodes (explicit nodes) comprise the initial inactive node set 
            do inode = 1, nnodes_total

                if (.not. (ANY(base_node_indices==inode))) then

                    call inactive_nodes%push_back(inode)

                end if

            end do

            ! Adaptively add the explicit nodes to the RBFI
            ! by moving one node at a time from the inactive to active node set.
            ! The inactive node with the largest distance to the active node set
            ! is selected to be added to the RBFI and moved to the active node set.

            num_inactive = inactive_nodes%size()

            num_active = active_nodes%size()

            do while (num_inactive>0)

                ! Loop over inactive nodes and find their distance to the active node set.
                ! Move the point with the greatest distance to the active node set.  
                
                dist_max = ZERO 
                inode_transfer = -1
                do inode = 1, num_inactive

                    inode_ig = inactive_nodes%at(inode)

                    
                    ! Find the distance to the closest active node
                    min_dist = 1.0e16_rk
                    do inode_active = 1, num_active


                        inode_ag = active_nodes%at(inode_active)

                        current_dist = node_dist(self%rbf_node_patch%nodes(inode_ig,:),self%rbf_node_patch%nodes(inode_ag,:))
                        
                        if (current_dist<min_dist) then
                            min_dist = current_dist
                        end if


                    end do

                    ! If this distance is the largest, mark the current inactive node to be activated
                    if (min_dist>dist_max) then
                        dist_max = min_dist
                        inode_transfer = inode_ig
                    end if

                end do
                radius_new = dist_max
                call rbf_node_temp%set_node(self%rbf_node_patch%nodes(inode_transfer,:), radius_new)
                call self%source_nodes_explicit%push_back(rbf_node_temp)
                call self%source_nodes%push_back(rbf_node_temp)

                call active_nodes%push_back(inode_transfer)
                call inactive_nodes%remove(inode_transfer)


                num_active = active_nodes%size()
                num_inactive = inactive_nodes%size()
                self%rbf_node_patch%patch_to_rbf_index(inode_transfer) = num_active
                self%rbf_node_patch%rbf_to_patch_index(num_active) = inode_transfer
                
            end do

            !
            ! Now, construct the interpolation matrix
            !

            call self%interpolation_matrix%init(nnodes_base, nnodes_explicit)


            do inode_source = 1, nnodes_base

                source_node = self%source_nodes_base%data(inode_source)%node_center
                radius_rbf = self%source_nodes_base%data(inode_source)%node_radius
                do inode_target = 1, nnodes_base

                    target_node = self%source_nodes_base%data(inode_target)%node_center
                    self%interpolation_matrix%A(inode_target, inode_source) = self%rbf%compute(target_node, source_node, radius_rbf)

                end do
                
                do inode_target = 1, nnodes_explicit

                    target_node = self%source_nodes_explicit%data(inode_target)%node_center
                    self%interpolation_matrix%B(inode_target, inode_source) = self%rbf%compute(target_node, source_node, radius_rbf)

                end do

            end do

            self%interpolation_matrix%Ainv = inv(self%interpolation_matrix%A)
            do inode_target = 1, nnodes_explicit

                target_node = self%source_nodes_explicit%data(inode_target)%node_center
                do inode_source = 1, inode_target
                    source_node = self%source_nodes_explicit%data(inode_source)%node_center
                    radius_rbf = self%source_nodes_explicit%data(inode_source)%node_radius

                    if (allocated(self%rbf_explicit)) then
                        self%interpolation_matrix%C(inode_target, inode_source) = self%rbf_explicit%compute(target_node, source_node, radius_rbf)
                    else
                        self%interpolation_matrix%C(inode_target, inode_source) = self%rbf%compute(target_node, source_node, radius_rbf)
                    end if

                end do
                
            end do
 
            self%interpolation_matrix%is_filled = .true.
        end if

    end subroutine construct_rbf_interpolation
    !********************************************************************************

    subroutine construct_rbf_interpolation_std(self, radius_base)
        class(rbf_interpolation_t),         intent(inout)   :: self
        real(rk),                           intent(in)      :: radius_base(3)

        type(ivector_t) :: active_nodes, inactive_nodes

        real(rk) :: radius_new(3), target_node(3), source_node(3), radius_rbf(3)
        real(rk) :: dist_max, min_dist, current_dist
        integer(ik) :: nnodes_total, nnodes_explicit
        integer(ik) :: inode, inode_active, inode_transfer, inode_ag, inode_ig, num_active, num_inactive, inode_source, inode_target

        type(rbf_node_t) :: rbf_node_temp

        integer(ik), allocatable :: base_node_indices(:)

        ! We may want to implement a NON-multiscale version to have a basis for comparision.

            ! Select the nnodes_base number of base nodes using k-means clustering
            nnodes_total = size(self%rbf_node_patch%nodes, 1)

            
            ! Store the base node informationinto the relevant data strucutres
            ! These comprise the initial active node set
            do inode = 1, nnodes_total

                call rbf_node_temp%set_node(self%rbf_node_patch%nodes(inode,:), radius_base)
                call active_nodes%push_back(inode)
                call self%source_nodes_base%push_back(rbf_node_temp)
                call self%source_nodes%push_back(rbf_node_temp)
                self%rbf_node_patch%patch_to_rbf_index(inode) = inode
                self%rbf_node_patch%rbf_to_patch_index(inode) = inode

            end do

            !
            ! Now, construct the interpolation matrix
            !

            allocate(self%interpolation_matrix%A(nnodes_total, nnodes_total))


            do inode_source = 1, nnodes_total

                source_node = self%source_nodes_base%data(inode_source)%node_center
                radius_rbf = self%source_nodes_base%data(inode_source)%node_radius
                do inode_target = 1, nnodes_total

                    target_node = self%source_nodes_base%data(inode_target)%node_center
                    self%interpolation_matrix%A(inode_target, inode_source) = self%rbf%compute(target_node, source_node, radius_rbf)

                end do
            end do

            self%interpolation_matrix%Ainv = inv(self%interpolation_matrix%A)

            self%interpolation_matrix%is_filled = .true.

    end subroutine construct_rbf_interpolation_std
    !********************************************************************************


    !>  Solves the RBF system with given source node values.
    !!
    !!  NOTE: The node ordering of source_node_vals is assumed to be the same as
    !!          the node ordering of rbf_node_patch%nodes. This is convenient, since
    !!          the we can easily evaluate PMM fomulas at the RBF patch nodes.
    !!  
    !!        The rbf_node_patch handles the reordering of these nodes to be
    !!          consistent with the RBF system node ordering.
    !! 
    !!  @author Eric Wolf
    !!  @date 10/19/2017
    !--------------------------------------------------------------------------------
    function solve(self, source_node_vals) result(rbf_coefficients)
        class(rbf_interpolation_t),         intent(inout)   :: self
        real(rk),                           intent(in)      :: source_node_vals(:)

        real(rk), allocatable   :: rbf_coefficients(:)

        real(rk), allocatable :: source_node_vals_rbf(:)

        !Re-index the source node values to conform to the RBF node ordering
        source_node_vals_rbf = self%rbf_node_patch%reorder_patch_to_rbf(source_node_vals)

        ! Now solve
        rbf_coefficients = self%interpolation_matrix%solve(source_node_vals_rbf) 


    end function solve
    !********************************************************************************

    !>  Evaluates a single given RBF at a given node.
    !!  
    !! 
    !!  @author Eric Wolf
    !!  @date 10/19/2017
    !--------------------------------------------------------------------------------
    function evaluate_rbf(self, inode, eval_node) result(val)
        class(rbf_interpolation_t),         intent(inout)  :: self
        integer(ik),                        intent(in)  :: inode
        real(rk),                           intent(in)  :: eval_node(3)

        real(rk)                                        :: val
        val = self%rbf%compute(eval_node, self%source_nodes%data(inode)%node_center, self%source_nodes%data(inode)%node_radius)

    end function evaluate_rbf
    !********************************************************************************

    !>  evaluates the gradient of a single given RBF at a given node.
    !!  
    !! 
    !!  @author Eric Wolf
    !!  @date 10/19/2017
    !--------------------------------------------------------------------------------
    function evaluate_rbf_grad(self, inode, eval_node) result(val)
        class(rbf_interpolation_t),         intent(inout)  :: self
        integer(ik),                        intent(in)  :: inode
        real(rk),                           intent(in)  :: eval_node(3)
        real(rk)                                        :: val(3)

        val = self%rbf%compute_grad(eval_node, self%source_nodes%data(inode)%node_center, self%source_nodes%data(inode)%node_radius)

    end function evaluate_rbf_grad
    !********************************************************************************

    !>  Evaluate the RBF interpolation at a given node.
    !!  
    !! 
    !!  @author Eric Wolf
    !!  @date 10/19/2017
    !--------------------------------------------------------------------------------
    function evaluate(self,rbf_coeff, eval_node) result(val)
        class(rbf_interpolation_t),         intent(inout)  :: self
        real(rk),                           intent(in)  :: rbf_coeff(:)
        real(rk),                           intent(in)  :: eval_node(3)

        real(rk)                                        :: val

        integer(ik)     :: isource_node

        ! Naive implementation, not using compact support
        val = ZERO
        do isource_node = 1, self%rbf_node_patch%nnodes_patch

            val = val + rbf_coeff(isource_node)*self%evaluate_rbf(isource_node,eval_node)

        end do
    end function evaluate
    !********************************************************************************

    !>  Evaluates the gradient of the RBF interpolant at a given node. 
    !!  
    !! 
    !!  @author Eric Wolf
    !!  @date 10/19/2017
    !--------------------------------------------------------------------------------
    function evaluate_grad(self,rbf_coeff, eval_node) result(val)
        class(rbf_interpolation_t),         intent(inout)  :: self
        real(rk),                           intent(in)  :: rbf_coeff(:)
        real(rk),                           intent(in)  :: eval_node(3)

        real(rk)                                        :: val(3)

        integer(ik)     :: isource_node

        ! Naive implementation, not using compact support

        val = ZERO
        do isource_node = 1, self%source_nodes%size()

            val = val + rbf_coeff(isource_node)*self%evaluate_rbf_grad(isource_node,eval_node)

        end do
    end function evaluate_grad
    !********************************************************************************

    !>  Evaluate the RBF interpolation at a given node.
    !!  
    !! 
    !!  @author Eric Wolf
    !!  @date 10/19/2017
    !--------------------------------------------------------------------------------
    function evaluate_ms(self,rbf_coeff, eval_node) result(val)
        class(rbf_interpolation_t),         intent(inout)  :: self
        real(rk),                           intent(in)  :: rbf_coeff(:)
        real(rk),                           intent(in)  :: eval_node(3)

        real(rk)                                        :: val

        integer(ik)     :: isource_node

        ! Naive implementation, not using compact support

        if (allocated(self%rbf_explicit)) then
            val = ZERO
            do isource_node = 1, self%source_nodes_base%size()

                val = val + rbf_coeff(isource_node)*self%rbf%compute(eval_node, &
                    self%source_nodes_base%data(isource_node)%node_center,&
                    self%source_nodes_base%data(isource_node)%node_radius)

            end do
            do isource_node = 1, self%source_nodes_explicit%size()

                val = val + rbf_coeff(self%source_nodes_base%size()+isource_node)*self%rbf_explicit%compute(eval_node, &
                    self%source_nodes_explicit%data(isource_node)%node_center,&
                    self%source_nodes_explicit%data(isource_node)%node_radius)

            end do

        else

            val = ZERO
            do isource_node = 1, self%rbf_node_patch%nnodes_patch

                val = val + rbf_coeff(isource_node)*self%evaluate_rbf(isource_node,eval_node)

            end do

        end if
        
    end function evaluate_ms
    !********************************************************************************

    !>  Evaluates the gradient of the RBF interpolant at a given node. 
    !!  
    !! 
    !!  @author Eric Wolf
    !!  @date 10/19/2017
    !--------------------------------------------------------------------------------
    function evaluate_ms_grad(self,rbf_coeff, eval_node) result(val)
        class(rbf_interpolation_t),         intent(inout)  :: self
        real(rk),                           intent(in)  :: rbf_coeff(:)
        real(rk),                           intent(in)  :: eval_node(3)

        real(rk)                                        :: val(3)

        integer(ik)     :: isource_node

        ! Naive implementation, not using compact support

        if (allocated(self%rbf_explicit)) then
            val = ZERO
            do isource_node = 1, self%source_nodes_base%size()

                val = val + rbf_coeff(isource_node)*self%rbf%compute_grad(eval_node, &
                    self%source_nodes_base%data(isource_node)%node_center,&
                    self%source_nodes_base%data(isource_node)%node_radius)

            end do
            do isource_node = 1, self%source_nodes_explicit%size()

                val = val + rbf_coeff(self%source_nodes_base%size()+isource_node)*self%rbf_explicit%compute_grad(eval_node, &
                    self%source_nodes_explicit%data(isource_node)%node_center,&
                    self%source_nodes_explicit%data(isource_node)%node_radius)

            end do

        else
            val = ZERO
            do isource_node = 1, self%source_nodes%size()

                val = val + rbf_coeff(isource_node)*self%evaluate_rbf_grad(isource_node,eval_node)

            end do

        end if

    end function evaluate_ms_grad
    !********************************************************************************



end module type_rbf_interpolation
