module mod_nodes_uniform
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: ZERO, ONE, TWO, THREE, FOUR
    use mod_gridspace,  only: linspace
    implicit none




contains




    !>
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/28/2017
    !!
    !-------------------------------------------------------------------------------
    function uniform_nodes(level, dim) result(nodes_)
        integer(ik),    intent(in)  :: level
        integer(ik),    intent(in)  :: dim

        ! Level for quadrature nodes corresponds to the number of nodes
        ! in the 1D tensor construction of the node set.
        real(rk),   dimension(level+1)          :: xi_nodes, eta_nodes, zeta_nodes
        integer(ik)                             :: ixi, ieta, izeta, inode, nnodes, ierr, nnodes1d
        integer(ik)                             :: nnodes_xi, nnodes_eta, nnodes_zeta
        real(rk),   allocatable, dimension(:)   :: nodes_(:,:)


        !
        ! Uniform distribution starts with two nodes at end points.
        ! Therefore, minimum level(1) is two points(level+1)
        !
        nnodes1d = level + 1


        !
        ! Get 1D nodes
        !
        select case(dim)
            case(1)
                nnodes_xi   = nnodes1d
                nnodes_eta  = 1
                nnodes_zeta = 1
                xi_nodes    = linspace(-ONE,ONE,nnodes1d)
                eta_nodes   = ZERO
                zeta_nodes  = ZERO
            case(2)
                nnodes_xi   = nnodes1d
                nnodes_eta  = nnodes1d
                nnodes_zeta = 1
                xi_nodes    = linspace(-ONE,ONE,nnodes1d)
                eta_nodes   = linspace(-ONE,ONE,nnodes1d)
                zeta_nodes  = ZERO
            case(3)
                nnodes_xi   = nnodes1d
                nnodes_eta  = nnodes1d
                nnodes_zeta = nnodes1d
                xi_nodes    = linspace(-ONE,ONE,nnodes1d)
                eta_nodes   = linspace(-ONE,ONE,nnodes1d)
                zeta_nodes  = linspace(-ONE,ONE,nnodes1d)
            case default
                call chidg_signal_one(FATAL,"uniform_nodes: invalid dimension.",dim)
        end select


        !
        ! Allocate nodes_
        !
        nnodes = nnodes_xi*nnodes_eta*nnodes_zeta
        allocate(nodes_(nnodes,3), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Contruct node set 
        !
        inode = 1
        do izeta = 1,nnodes_zeta
            do ieta = 1,nnodes_eta
                do ixi = 1,nnodes_xi
                    nodes_(inode,1) = xi_nodes(ixi)
                    nodes_(inode,2) = eta_nodes(ieta)
                    nodes_(inode,3) = zeta_nodes(izeta)
                    inode = inode + 1
                end do
            end do
        end do


    end function uniform_nodes
    !***********************************************************************************







    !>
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/5/2017
    !!
    !-----------------------------------------------------------------------------------
    function uniform_weights(level,dim) result(weights)
        integer(ik),    intent(in)  :: level
        integer(ik),    intent(in)  :: dim

        real(rk),   allocatable :: weights(:), xi_weights(:), eta_weights(:), zeta_weights(:)
        integer(ik)             :: nnodes1d, nweights_xi, nweights_eta, nweights_zeta, nweights, &
                                   inode, ixi, ieta, izeta, ierr
        
        !
        ! For the node set 'Uniform', the number of nodes in 1D is
        ! defined as (level + 1)
        !
        nnodes1d = level + 1



        !
        ! Get 1D weights
        !
        allocate(xi_weights(nnodes1d), eta_weights(nnodes1d), zeta_weights(nnodes1d), stat=ierr)
        if (ierr /= 0) call AllocationError

        select case(dim)
            case(1)
                nweights_xi   = nnodes1d
                nweights_eta  = 1
                nweights_zeta = 1
                xi_weights    = uniform_weights_1d(nnodes1d)
                eta_weights   = ONE
                zeta_weights  = ONE
            case(2)
                nweights_xi   = nnodes1d
                nweights_eta  = nnodes1d
                nweights_zeta = 1
                xi_weights    = uniform_weights_1d(nnodes1d)
                eta_weights   = uniform_weights_1d(nnodes1d)
                zeta_weights  = ONE
            case(3)
                nweights_xi   = nnodes1d
                nweights_eta  = nnodes1d
                nweights_zeta = nnodes1d
                xi_weights    = uniform_weights_1d(nnodes1d)
                eta_weights   = uniform_weights_1d(nnodes1d)
                zeta_weights  = uniform_weights_1d(nnodes1d)
            case default
                call chidg_signal_one(FATAL,"quadrature_weights: invalid dimension.",dim)
        end select


        !
        ! Allocate weights_
        !
        nweights = nweights_xi*nweights_eta*nweights_zeta
        allocate(weights(nweights), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Contruct node set 
        !
        inode = 1
        do izeta = 1,nweights_zeta
            do ieta = 1,nweights_eta
                do ixi = 1,nweights_xi
                    weights(inode) = xi_weights(ixi)*eta_weights(ieta)*zeta_weights(izeta)
                    inode = inode + 1
                end do
            end do
        end do

    end function uniform_weights
    !***********************************************************************************




















!    !>
!    !!
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   6/28/2017
!    !!
!    !-------------------------------------------------------------------------------
!    function uniform_nodes(level, dim) result(nodes)
!        integer(ik),    intent(in)  :: level
!        integer(ik),    intent(in)  :: dim
!
!        ! Level for quadrature nodes corresponds to the number of nodes
!        ! in the 1D tensor construction of the node set.
!        integer(ik)                 :: inode, inode_inner, nnodes, ierr, &
!                                       nnodes1d, nedge_vertices, nedge_interiors, nface_interiors, &
!                                       nelement_interiors, node_start, node_end, line_begin, line_end
!        real(rk),       allocatable :: nodes(:,:), line_space(:)
!        integer(ik),    allocatable :: nodes_1daccess(:,:), interior_forward(:), interior_reverse(:)
!        logical                     :: duplicate_node
!
!
!        if (dim /= 3) call chidg_signal_one(FATAL,'uniform_nodes: only dim=3 is valid.',dim)
!
!
!        nnodes1d = level+1
!        line_space = linspace(-ONE,ONE,nnodes1d)
!
!        !
!        ! Uniform distribution starts with two nodes at end points.
!        ! Therefore, minimum level(1) is two points(level+1)
!        !
!        nnodes1d = level + 1
!        nnodes = nnodes1d*nnodes1d*nnodes1d
!        allocate(nodes(nnodes,3), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!
!        ! Get tensor product access indices
!        nodes_1daccess = get_reference_element_1daccess(level)
!
!        ! Check dimensions are correct
!        if (size(nodes_1daccess,1) /= size(nodes,1)) call chidg_signal(FATAL,'uniform_weights: inconsistent node dimensions.')
!
!        ! Assemble nodes
!        do inode = 1,size(nodes,1)
!            nodes(inode,1) = line_space(nodes_1daccess(inode,1))
!            nodes(inode,2) = line_space(nodes_1daccess(inode,2))
!            nodes(inode,3) = line_space(nodes_1daccess(inode,3))
!        end do !inode
!
!
!        ! Check to make sure there are no duplicate nodes
!        do inode = 1,size(nodes,1)
!
!            duplicate_node = .false.
!            do inode_inner = 1,size(nodes,1)
!                if ( (inode /= inode_inner) .and. &
!                     (nodes(inode,1) == nodes(inode_inner,1)) .and. &
!                     (nodes(inode,2) == nodes(inode_inner,2)) .and. &
!                     (nodes(inode,3) == nodes(inode_inner,3)) &
!                    ) then
!                    duplicate_node = .true. 
!                end if
!            end do !inode_inner
!
!            if (duplicate_node) call chidg_signal(FATAL,'uniform_nodes: found duplicate node. Must be error in node set convention.')
!
!        end do !inode
!
!
!
!
!    end function uniform_nodes
!    !***********************************************************************************
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!    !>
!    !!
!    !!
!    !!  @author Nathan A. Wukie (AFRL)
!    !!  @date   7/5/2017
!    !!
!    !-----------------------------------------------------------------------------------
!    function uniform_weights(level,dim) result(weights)
!        integer(ik),    intent(in)  :: level
!        integer(ik),    intent(in)  :: dim
!
!        real(rk),       allocatable :: weights(:), line_weights(:)
!        integer(ik),    allocatable :: nodes_1daccess(:,:)
!        integer(ik)                 :: nnodes1d, nweights_xi, nweights_eta, nweights_zeta, nweights, &
!                                       inode, ierr, iweight
!        
!        ! For the node set 'Uniform', the number of nodes in 1D is
!        ! defined as (level + 1)
!        nnodes1d = level + 1
!
!
!        ! Get line weights
!        line_weights = uniform_weights_1d(nnodes1d)
!
!        ! Get tensor product access indices
!        nodes_1daccess = get_reference_element_1daccess(level)
!
!        ! Specialize for 
!        select case(dim)
!            case(3)
!                nweights_xi   = nnodes1d
!                nweights_eta  = nnodes1d
!                nweights_zeta = nnodes1d
!            case default
!                call chidg_signal_one(FATAL,"quadrature_weights: invalid dimension.",dim)
!        end select
!
!        ! Allocate weights_
!        nweights = nweights_xi*nweights_eta*nweights_zeta
!        allocate(weights(nweights), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!
!        ! Check dimensions are correct
!        if (size(nodes_1daccess,1) /= size(weights)) call chidg_signal(FATAL,'uniform_weights: inconsistent node dimensions.')
!
!        
!        ! Assemble tensor product weights
!        do iweight = 1,size(weights)
!            weights(iweight) = line_weights(nodes_1daccess(iweight,1)) * &
!                               line_weights(nodes_1daccess(iweight,2)) * &
!                               line_weights(nodes_1daccess(iweight,3))
!        end do !iweight
!
!    end function uniform_weights
!    !***********************************************************************************
!
!
!
!
!
!
!
!
!
!
!    !>
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   8/1/2019
!    !!
!    !----------------------------------------------------------------------------------
!    function get_reference_element_1daccess(order) result(nodes_1daccess)
!        integer(ik),    intent(in)  :: order
!
!        integer(ik)                 :: nedge_vertices, nedge_interiors, nface_interiors, nelement_interiors, &
!                                       line_begin, line_end, node_start, node_end, nnodes, nnodes1d, ierr, inode, inode_inner
!        integer(ik),    allocatable :: nodes_1daccess(:,:), interior_forward(:), interior_reverse(:)
!        logical :: duplicate_node
!
!
!        nnodes1d = order+1
!        nnodes = nnodes1d*nnodes1d*nnodes1d
!        allocate(nodes_1daccess(nnodes,3), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!
!        !
!        ! Fill edge, face, interior nodes
!        !
!        select case (order)
!            ! Quadratic (HEXA_27)
!            case(2)
!                nedge_vertices     = 2
!                nedge_interiors    = 1
!                nface_interiors    = 1
!                nelement_interiors = 1
!                line_begin = 1
!                line_end  = nedge_interiors+nedge_vertices
!                interior_forward = [2]
!                interior_reverse = [2]
!
!
!
!
!                ! Hex vertices are same for all orders
!                nodes_1daccess(1,1) = line_begin
!                nodes_1daccess(1,2) = line_begin
!                nodes_1daccess(1,3) = line_begin
!
!                nodes_1daccess(2,1) = line_end
!                nodes_1daccess(2,2) = line_begin
!                nodes_1daccess(2,3) = line_begin
!
!                nodes_1daccess(3,1) = line_end
!                nodes_1daccess(3,2) = line_end
!                nodes_1daccess(3,3) = line_begin
!
!                nodes_1daccess(4,1) = line_begin
!                nodes_1daccess(4,2) = line_end
!                nodes_1daccess(4,3) = line_begin
!
!                nodes_1daccess(5,1) = line_begin
!                nodes_1daccess(5,2) = line_begin
!                nodes_1daccess(5,3) = line_end
!
!                nodes_1daccess(6,1) = line_end
!                nodes_1daccess(6,2) = line_begin
!                nodes_1daccess(6,3) = line_end
!
!                nodes_1daccess(7,1) = line_end
!                nodes_1daccess(7,2) = line_end
!                nodes_1daccess(7,3) = line_end
!
!                nodes_1daccess(8,1) = line_begin
!                nodes_1daccess(8,2) = line_end
!                nodes_1daccess(8,3) = line_end
!
!
!
!
!
!
!
!
!
!
!
!                ! Edge interiors
!
!                !* Edge 1
!                node_start = 9
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = interior_forward
!                nodes_1daccess(node_start:node_end, 2) = line_begin
!                nodes_1daccess(node_start:node_end, 3) = line_begin
!
!                !* Edge 2
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_end
!                nodes_1daccess(node_start:node_end, 2) = interior_forward
!                nodes_1daccess(node_start:node_end, 3) = line_begin
!
!                !* Edge 3
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = interior_reverse
!                nodes_1daccess(node_start:node_end, 2) = line_end
!                nodes_1daccess(node_start:node_end, 3) = line_begin
!
!                !* Edge 4
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_begin
!                nodes_1daccess(node_start:node_end, 2) = interior_reverse
!                nodes_1daccess(node_start:node_end, 3) = line_begin
!
!                !* Edge 5
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_begin
!                nodes_1daccess(node_start:node_end, 2) = line_begin
!                nodes_1daccess(node_start:node_end, 3) = interior_forward
!
!                !* Edge 6
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_end
!                nodes_1daccess(node_start:node_end, 2) = line_begin
!                nodes_1daccess(node_start:node_end, 3) = interior_forward
!
!                !* Edge 7
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_end
!                nodes_1daccess(node_start:node_end, 2) = line_end
!                nodes_1daccess(node_start:node_end, 3) = interior_forward
!
!                !* Edge 8
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_begin
!                nodes_1daccess(node_start:node_end, 2) = line_end
!                nodes_1daccess(node_start:node_end, 3) = interior_forward
!
!                !* Edge 9
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = interior_forward
!                nodes_1daccess(node_start:node_end, 2) = line_begin
!                nodes_1daccess(node_start:node_end, 3) = line_end
!
!                !* Edge 10
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_end
!                nodes_1daccess(node_start:node_end, 2) = interior_forward
!                nodes_1daccess(node_start:node_end, 3) = line_end
!
!                !* Edge 11
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = interior_reverse
!                nodes_1daccess(node_start:node_end, 2) = line_end
!                nodes_1daccess(node_start:node_end, 3) = line_end
!
!                !* Edge 12
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_begin
!                nodes_1daccess(node_start:node_end, 2) = interior_reverse
!                nodes_1daccess(node_start:node_end, 3) = line_end
!
!
!
!                ! Face interiors
!
!                ! Face 1
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end,1) = 2
!                nodes_1daccess(node_start:node_end,2) = 2
!                nodes_1daccess(node_start:node_end,3) = line_begin
!
!                ! Face 2
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end,1) = 2
!                nodes_1daccess(node_start:node_end,2) = line_begin
!                nodes_1daccess(node_start:node_end,3) = 2
!
!                ! Face 3
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end,1) = line_end
!                nodes_1daccess(node_start:node_end,2) = 2
!                nodes_1daccess(node_start:node_end,3) = 2
!
!                ! Face 4
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end,1) = 2
!                nodes_1daccess(node_start:node_end,2) = line_end
!                nodes_1daccess(node_start:node_end,3) = 2
!
!                ! Face 5
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end,1) = line_begin
!                nodes_1daccess(node_start:node_end,2) = 2
!                nodes_1daccess(node_start:node_end,3) = 2
!
!                ! Face 6
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end,1) = 2
!                nodes_1daccess(node_start:node_end,2) = 2
!                nodes_1daccess(node_start:node_end,3) = line_end
!
!                !* Interior
!                node_start = node_end+1
!                node_end   = node_start+nelement_interiors-1
!                nodes_1daccess(node_start:node_end,1) = 2
!                nodes_1daccess(node_start:node_end,2) = 2
!                nodes_1daccess(node_start:node_end,3) = 2
!
!
!
!            ! Cubic (HEXA_64)
!            case(3)
!
!                nedge_vertices     = 2
!                nedge_interiors    = 2
!                nface_interiors    = 4
!                nelement_interiors = 8
!                line_begin = 1
!                line_end  = nedge_interiors+nedge_vertices
!                interior_forward = [2,3]
!                interior_reverse = [3,2]
!
!
!
!                ! Hex vertices are same for all orders
!                nodes_1daccess(1,1) = line_begin
!                nodes_1daccess(1,2) = line_begin
!                nodes_1daccess(1,3) = line_begin
!
!                nodes_1daccess(2,1) = line_end
!                nodes_1daccess(2,2) = line_begin
!                nodes_1daccess(2,3) = line_begin
!
!                nodes_1daccess(3,1) = line_end
!                nodes_1daccess(3,2) = line_end
!                nodes_1daccess(3,3) = line_begin
!
!                nodes_1daccess(4,1) = line_begin
!                nodes_1daccess(4,2) = line_end
!                nodes_1daccess(4,3) = line_begin
!
!                nodes_1daccess(5,1) = line_begin
!                nodes_1daccess(5,2) = line_begin
!                nodes_1daccess(5,3) = line_end
!
!                nodes_1daccess(6,1) = line_end
!                nodes_1daccess(6,2) = line_begin
!                nodes_1daccess(6,3) = line_end
!
!                nodes_1daccess(7,1) = line_end
!                nodes_1daccess(7,2) = line_end
!                nodes_1daccess(7,3) = line_end
!
!                nodes_1daccess(8,1) = line_begin
!                nodes_1daccess(8,2) = line_end
!                nodes_1daccess(8,3) = line_end
!
!
!
!
!
!
!
!
!
!                ! Edge interiors
!                !* Edge 1
!                node_start = 9
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = interior_forward
!                nodes_1daccess(node_start:node_end, 2) = line_begin
!                nodes_1daccess(node_start:node_end, 3) = line_begin
!
!                !* Edge 2
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_end
!                nodes_1daccess(node_start:node_end, 2) = interior_forward
!                nodes_1daccess(node_start:node_end, 3) = line_begin
!
!                !* Edge 3
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = interior_reverse
!                nodes_1daccess(node_start:node_end, 2) = line_end
!                nodes_1daccess(node_start:node_end, 3) = line_begin
!
!                !* Edge 4
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_begin
!                nodes_1daccess(node_start:node_end, 2) = interior_reverse
!                nodes_1daccess(node_start:node_end, 3) = line_begin
!
!                !* Edge 5
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_begin
!                nodes_1daccess(node_start:node_end, 2) = line_begin
!                nodes_1daccess(node_start:node_end, 3) = interior_forward
!
!                !* Edge 6
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_end
!                nodes_1daccess(node_start:node_end, 2) = line_begin
!                nodes_1daccess(node_start:node_end, 3) = interior_forward
!
!                !* Edge 7
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_end
!                nodes_1daccess(node_start:node_end, 2) = line_end
!                nodes_1daccess(node_start:node_end, 3) = interior_forward
!
!                !* Edge 8
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_begin
!                nodes_1daccess(node_start:node_end, 2) = line_end
!                nodes_1daccess(node_start:node_end, 3) = interior_forward
!
!                !* Edge 9
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = interior_forward
!                nodes_1daccess(node_start:node_end, 2) = line_begin
!                nodes_1daccess(node_start:node_end, 3) = line_end
!
!                !* Edge 10
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_end
!                nodes_1daccess(node_start:node_end, 2) = interior_forward
!                nodes_1daccess(node_start:node_end, 3) = line_end
!
!                !* Edge 11
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = interior_reverse
!                nodes_1daccess(node_start:node_end, 2) = line_end
!                nodes_1daccess(node_start:node_end, 3) = line_end
!
!                !* Edge 12
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_begin
!                nodes_1daccess(node_start:node_end, 2) = interior_reverse
!                nodes_1daccess(node_start:node_end, 3) = line_end
!
!
!
!                !* Face 1
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = [2, 3, 3, 2]
!                nodes_1daccess(node_start:node_end, 2) = [2, 2, 3, 3]
!                nodes_1daccess(node_start:node_end, 3) = line_begin
!
!                !* Face 2
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = [2, 3, 3, 2]
!                nodes_1daccess(node_start:node_end, 2) = line_begin
!                nodes_1daccess(node_start:node_end, 3) = [2, 2, 3, 3]
!
!                !* Face 3
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_end
!                nodes_1daccess(node_start:node_end, 2) = [2, 3, 3, 2]
!                nodes_1daccess(node_start:node_end, 3) = [2, 2, 3, 3]
!
!                !* Face 4
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = [3, 2, 2, 3]
!                nodes_1daccess(node_start:node_end, 2) = line_end
!                nodes_1daccess(node_start:node_end, 3) = [2, 2, 3, 3]
!
!                !* Face 5
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_begin
!                nodes_1daccess(node_start:node_end, 2) = [3, 2, 2, 3]
!                nodes_1daccess(node_start:node_end, 3) = [2, 2, 3, 3]
!
!                !* Face 6
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = [2, 3, 3, 2]
!                nodes_1daccess(node_start:node_end, 2) = [2, 2, 3, 3]
!                nodes_1daccess(node_start:node_end, 3) = line_end
!
!                !* Interior
!                node_start = node_end+1
!                node_end   = node_start+nelement_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = [2, 3, 3, 2, 2, 3, 3, 2]
!                nodes_1daccess(node_start:node_end, 2) = [2, 2, 3, 3, 2, 2, 3, 3]
!                nodes_1daccess(node_start:node_end, 3) = [2, 2, 2, 2, 3, 3, 3, 3]
!
!
!
!            ! Quartic (HEXA_125)
!            case(4)
!
!                nedge_vertices     = 2
!                nedge_interiors    = 3
!                nface_interiors    = 9
!                nelement_interiors = 27
!                line_begin = 1
!                line_end  = nedge_interiors+nedge_vertices
!                interior_forward = [2,3,4]
!                interior_reverse = [4,3,2]
!
!
!
!                ! Hex vertices are same for all orders
!                nodes_1daccess(1,1) = line_begin
!                nodes_1daccess(1,2) = line_begin
!                nodes_1daccess(1,3) = line_begin
!
!                nodes_1daccess(2,1) = line_end
!                nodes_1daccess(2,2) = line_begin
!                nodes_1daccess(2,3) = line_begin
!
!                nodes_1daccess(3,1) = line_end
!                nodes_1daccess(3,2) = line_end
!                nodes_1daccess(3,3) = line_begin
!
!                nodes_1daccess(4,1) = line_begin
!                nodes_1daccess(4,2) = line_end
!                nodes_1daccess(4,3) = line_begin
!
!                nodes_1daccess(5,1) = line_begin
!                nodes_1daccess(5,2) = line_begin
!                nodes_1daccess(5,3) = line_end
!
!                nodes_1daccess(6,1) = line_end
!                nodes_1daccess(6,2) = line_begin
!                nodes_1daccess(6,3) = line_end
!
!                nodes_1daccess(7,1) = line_end
!                nodes_1daccess(7,2) = line_end
!                nodes_1daccess(7,3) = line_end
!
!                nodes_1daccess(8,1) = line_begin
!                nodes_1daccess(8,2) = line_end
!                nodes_1daccess(8,3) = line_end
!
!
!
!
!
!
!
!
!
!
!
!                ! Edge interiors
!                !* Edge 1
!                node_start = 9
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = interior_forward
!                nodes_1daccess(node_start:node_end, 2) = line_begin
!                nodes_1daccess(node_start:node_end, 3) = line_begin
!
!                !* Edge 2
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_end
!                nodes_1daccess(node_start:node_end, 2) = interior_forward
!                nodes_1daccess(node_start:node_end, 3) = line_begin
!
!                !* Edge 3
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = interior_reverse
!                nodes_1daccess(node_start:node_end, 2) = line_end
!                nodes_1daccess(node_start:node_end, 3) = line_begin
!
!                !* Edge 4
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_begin
!                nodes_1daccess(node_start:node_end, 2) = interior_reverse
!                nodes_1daccess(node_start:node_end, 3) = line_begin
!
!                !* Edge 5
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_begin
!                nodes_1daccess(node_start:node_end, 2) = line_begin
!                nodes_1daccess(node_start:node_end, 3) = interior_forward
!
!                !* Edge 6
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_end
!                nodes_1daccess(node_start:node_end, 2) = line_begin
!                nodes_1daccess(node_start:node_end, 3) = interior_forward
!
!                !* Edge 7
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_end
!                nodes_1daccess(node_start:node_end, 2) = line_end
!                nodes_1daccess(node_start:node_end, 3) = interior_forward
!
!                !* Edge 8
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_begin
!                nodes_1daccess(node_start:node_end, 2) = line_end
!                nodes_1daccess(node_start:node_end, 3) = interior_forward
!
!                !* Edge 9
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = interior_forward
!                nodes_1daccess(node_start:node_end, 2) = line_begin
!                nodes_1daccess(node_start:node_end, 3) = line_end
!
!                !* Edge 10
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_end
!                nodes_1daccess(node_start:node_end, 2) = interior_forward
!                nodes_1daccess(node_start:node_end, 3) = line_end
!
!                !* Edge 11
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = interior_reverse
!                nodes_1daccess(node_start:node_end, 2) = line_end
!                nodes_1daccess(node_start:node_end, 3) = line_end
!
!                !* Edge 12
!                node_start = node_end+1
!                node_end   = node_start+nedge_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_begin
!                nodes_1daccess(node_start:node_end, 2) = interior_reverse
!                nodes_1daccess(node_start:node_end, 3) = line_end
!
!
!
!                !* Face 1
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = [2, 3, 4, 4, 4, &
!                                                          3, 2, 2, 3]
!                nodes_1daccess(node_start:node_end, 2) = [2, 2, 2, 3, 4, &
!                                                          4, 4, 3, 3]
!                nodes_1daccess(node_start:node_end, 3) = line_begin
!
!
!                !* Face 2
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = [2, 3, 4, 4, 4, &
!                                                          3, 2, 2, 3]
!                nodes_1daccess(node_start:node_end, 2) = line_begin
!                nodes_1daccess(node_start:node_end, 3) = [2, 2, 2, 3, 4, &
!                                                          4, 4, 3, 3]
!
!                !* Face 3
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_end
!                nodes_1daccess(node_start:node_end, 2) = [2, 3, 4, 4, 4, &
!                                                          3, 2, 2, 3]
!                nodes_1daccess(node_start:node_end, 3) = [2, 2, 2, 3, 4, &
!                                                          4, 4, 3, 3]
!
!                !* Face 4
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = [4, 3, 2, 2, 2, &
!                                                          3, 4, 4, 3]
!                nodes_1daccess(node_start:node_end, 2) = line_end
!                nodes_1daccess(node_start:node_end, 3) = [2, 2, 2, 3, 4, &
!                                                          4, 4, 3, 3]
!
!                !* Face 5
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = line_begin
!                nodes_1daccess(node_start:node_end, 2) = [4, 3, 2, 2, 2, &
!                                                          3, 4, 4, 3]
!                nodes_1daccess(node_start:node_end, 3) = [2, 2, 2, 3, 4, &
!                                                          4, 4, 3, 3]
!
!                !* Face 6
!                node_start = node_end+1
!                node_end   = node_start+nface_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = [2, 3, 4, 4, 4, &
!                                                          3, 2, 2, 3]
!                nodes_1daccess(node_start:node_end, 2) = [2, 2, 2, 3, 4, &
!                                                          4, 4, 3, 3]
!                nodes_1daccess(node_start:node_end, 3) = line_end
!
!                !* Interior
!                node_start = node_end+1
!                node_end   = node_start+nelement_interiors-1
!                nodes_1daccess(node_start:node_end, 1) = [2, 3, 4, 4, 4, &
!                                                          3, 2, 2, 3,    &
!                                                          2, 3, 4, 4, 4, &
!                                                          3, 2, 2, 3,    &
!                                                          2, 3, 4, 4, 4, &
!                                                          3, 2, 2, 3]
!                nodes_1daccess(node_start:node_end, 2) = [2, 2, 2, 3, 4, &
!                                                          4, 4, 3, 3,    &
!                                                          2, 2, 2, 3, 4, &
!                                                          4, 4, 3, 3,    &
!                                                          2, 2, 2, 3, 4, &
!                                                          4, 4, 3, 3]
!                nodes_1daccess(node_start:node_end, 3) = [2, 2, 2, &
!                                                          2, 2, 2, & 
!                                                          2, 2, 2, & 
!                                                          3, 3, 3, & 
!                                                          3, 3, 3, & 
!                                                          3, 3, 3, & 
!                                                          4, 4, 4, & 
!                                                          4, 4, 4, & 
!                                                          4, 4, 4]
!
!
!
!
!            case default
!                call chidg_signal_one(FATAL,"uniform_nodes: invalid value for 'level'", order)
!
!        end select
!
!
!        ! Check access indices are all unique
!        do inode = 1,size(nodes_1daccess,1)
!            duplicate_node = .false.
!            do inode_inner = 1,size(nodes_1daccess,1)
!                if ( (inode /= inode_inner) .and. &
!                     (nodes_1daccess(inode,1) == nodes_1daccess(inode_inner,1)) .and. &
!                     (nodes_1daccess(inode,2) == nodes_1daccess(inode_inner,2)) .and. &
!                     (nodes_1daccess(inode,3) == nodes_1daccess(inode_inner,3)) ) then
!                     duplicate_node = .true.
!                end if 
!            end do
!            if (duplicate_node) call chidg_signal(FATAL,'get_reference_element_1daccess: duplicate node detected.')
!        end do
!
!
!
!    end function get_reference_element_1daccess
!    !**********************************************************************************

















































































    !>  Compute quadrature weights for a 1D uniform node distribution between [-1,1].
    !!
    !!  Even number of nodes: Composite Trapezoid rule (less accurate)
    !!  Odd number of nodes:  Composite Simpson's rule (more accurate)
    !!
    !!
    !!  NOTE: An odd number of nodes should always be preferred, since the Composite
    !!        Simpson's rule is much more accurate than the Composite Trapezoid rule.
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/5/2017
    !!
    !------------------------------------------------------------------------------------
    function uniform_weights_1d(nnodes1d) result(weights)
        integer(ik),    intent(in)  :: nnodes1d

        real(rk)    :: weights(nnodes1d)
        integer(ik) :: inode
        logical     :: even, even_node

        ! determin odd/even
        even = (mod(nnodes1d,2) == 0)


        !
        ! EVEN: Composite Trapezoid rule
        !   int_a^b[f]dx = h/2[f(x0) + 2f(x1) + 2f(x2) + ... + 2f(x_{n-1}) + f(xn)]
        !   Weights computed for x in [-1,1]. h=2
        !
        if (even) then
            weights(1)        = ONE ! f(x0)
            weights(nnodes1d) = ONE ! f(xn)

            
            if (nnodes1d > 2) then
                do inode = 2,nnodes1d-1
                    weights(inode) = TWO
                end do !inode
            end if


        !
        ! ODD: Composite Simpson's rule
        !   int_a^b[f]dx = h/3[f(x0) + 4f(x1) + 2f(x2) + 4f(x3) + 2f(x4) + ... + 4f(x_{n-1}) + f(xn)]
        !   Weights computed for x in [-1,1]. h=2
        !
        else
            weights(1)        = TWO/THREE   ! f(x0)
            weights(nnodes1d) = TWO/THREE   ! f(xn)

            do inode = 2,nnodes1d-1
                even_node = (mod(inode,2) == 0)
                if (even_node) then
                    weights(inode) = (TWO/THREE)*FOUR
                else
                    weights(inode) = (TWO/THREE)*TWO
                end if
            end do !inode

        end if


    end function uniform_weights_1d
    !************************************************************************************




end module mod_nodes_uniform
