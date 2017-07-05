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



















!    !>
!    !!
!    !!
!    !!  @author Nathan A. Wukie (AFRL)
!    !!  @date   6/28/2017
!    !!
!    !-----------------------------------------------------------------------------------
!    function linspace(xmin,xmax,res) result(nodes_)
!        real(rk),       intent(in)  :: xmin
!        real(rk),       intent(in)  :: xmax
!        integer(ik),    intent(in)  :: res
!        
!        integer(ik) :: inode
!        real(rk)    :: nodes_(res)
!
!        !
!        ! Compute 1d point coordinates in each direction
!        !
!        do inode = 1,res
!            nodes_(inode) = xmin + (xmax-xmin)*(real(inode,rk)-ONE)/(real(res,rk)-ONE)
!        end do
!
!    end function linspace
!    !***********************************************************************************

















end module mod_nodes_uniform
