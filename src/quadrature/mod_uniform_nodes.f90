module mod_nodes_uniform
    use mod_kinds,  only: rk, ik
    use mod_constants,  only: ZERO, ONE, TWO




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
        real(rk),   dimension(level)            :: xi_nodes, eta_nodes, zeta_nodes
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
                nnodes_xi   = nnode1d
                nnodes_eta  = 1
                nnodes_zeta = 1
                xi_nodes    = linspace(-ONE,ONE,nnodes1d)
                eta_nodes   = ZERO
                zeta_nodes  = ZERO
            case(2)
                nnodes_xi   = nnode1d
                nnodes_eta  = nnode1d
                nnodes_zeta = 1
                xi_nodes    = linspace(-ONE,ONE,nnodes1d)
                eta_nodes   = linspace(-ONE,ONE,nnodes1d)
                zeta_nodes  = ZERO
            case(3)
                nnodes_xi   = nnode1d
                nnodes_eta  = nnode1d
                nnodes_zeta = nnode1d
                xi_nodes    = linspace(-ONE,ONE,nnodes1d)
                eta_nodes   = linspace(-ONE,ONE,nnodes1d)
                zeta_nodes  = linspace(-ONE,ONE,nnodes1d)
            case default
                call chidg_signal_one(FATAL,"uniform_nodes: invalid dimension.",dim)
        end select

        print*, xi_nodes

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
    !!  @date   6/28/2017
    !!
    !-----------------------------------------------------------------------------------
    function linspace(xmin,xmax,res) result(nodes_)
        real(rk),       intent(in)  :: xmin
        real(rk),       intent(in)  :: xmax
        integer(ik),    intent(in)  :: res
        
        integer(ik)             :: inode
        real(rk),   allocatable :: nodes_(:)

        !
        ! Compute 1d point coordinates in each direction
        !
        do inode = 1,res
            nodes_(inode) = xmin + (xmax-xmin)*(real(inode,rk)-ONE)/(real(res,rk)-ONE)
        end do

    end function linspace
    !***********************************************************************************

















end module mod_nodes_uniform
