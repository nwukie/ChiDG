module mod_gauss_legendre
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: PI, ONE, TWO, THREE, FOUR, EIGHT, ZERO
    use mod_legendre,       only: LegendreVal1D, DLegendreVal1D
    implicit none

contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/28/2017
    !!
    !------------------------------------------------------------------------
    function quadrature_nodes(level,dim) result(nodes_)
        integer(ik),    intent(in)  :: level
        integer(ik),    intent(in)  :: dim

        ! Level for quadrature nodes corresponds to the number of nodes
        ! in the 1D tensor construction of the node set.
        real(rk),   dimension(level)            :: xi_nodes, eta_nodes, zeta_nodes
        integer(ik)                             :: ixi, ieta, izeta, inode, nnodes, ierr
        integer(ik)                             :: nnodes_xi, nnodes_eta, nnodes_zeta
        real(rk),   allocatable, dimension(:)   :: nodes_(:,:)

        !
        ! Get 1D nodes
        !
        select case(dim)
            case(1)
                nnodes_xi   = level
                nnodes_eta  = 1
                nnodes_zeta = 1
                call gl_nodes(level,xi_nodes)
                eta_nodes  = ZERO
                zeta_nodes = ZERO
            case(2)
                nnodes_xi   = level
                nnodes_eta  = level
                nnodes_zeta = 1
                call gl_nodes(level,xi_nodes)
                call gl_nodes(level,eta_nodes)
                zeta_nodes = ZERO
            case(3)
                nnodes_xi   = level
                nnodes_eta  = level
                nnodes_zeta = level
                call gl_nodes(level,xi_nodes)
                call gl_nodes(level,eta_nodes)
                call gl_nodes(level,zeta_nodes)
            case default
                call chidg_signal_one(FATAL,"quadrature_nodes: invalid dimension.",dim)
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


    end function quadrature_nodes
    !************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/28/2017
    !!
    !------------------------------------------------------------------------
    function quadrature_weights(level,dim) result(weights_)
        integer(ik),    intent(in)  :: level
        integer(ik),    intent(in)  :: dim

        ! Level for quadrature weights corresponds to the number of weights
        ! in the 1D tensor construction of the node set.
        real(rk),   dimension(level)            :: xi_weights, eta_weights, zeta_weights
        integer(ik)                             :: ixi, ieta, izeta, inode, nweights, ierr
        integer(ik)                             :: nweights_xi, nweights_eta, nweights_zeta
        real(rk),   allocatable, dimension(:)   :: weights_(:)


        !
        ! Get 1D weights
        !
        select case(dim)
            case(1)
                nweights_xi   = level
                nweights_eta  = 1
                nweights_zeta = 1
                call gl_weights(level,xi_weights)
                eta_weights  = ONE
                zeta_weights = ONE
            case(2)
                nweights_xi   = level
                nweights_eta  = level
                nweights_zeta = 1
                call gl_weights(level,xi_weights)
                call gl_weights(level,eta_weights)
                zeta_weights = ONE
            case(3)
                nweights_xi   = level
                nweights_eta  = level
                nweights_zeta = level
                call gl_weights(level,xi_weights)
                call gl_weights(level,eta_weights)
                call gl_weights(level,zeta_weights)
            case default
                call chidg_signal_one(FATAL,"quadrature_weights: invalid dimension.",dim)
        end select


        !
        ! Allocate weights_
        !
        nweights = nweights_xi*nweights_eta*nweights_zeta
        allocate(weights_(nweights), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Contruct node set 
        !
        inode = 1
        do izeta = 1,nweights_zeta
            do ieta = 1,nweights_eta
                do ixi = 1,nweights_xi
                    weights_(inode) = xi_weights(ixi)*eta_weights(ieta)*zeta_weights(izeta)
                    inode = inode + 1
                end do
            end do
        end do


    end function quadrature_weights
    !************************************************************************









    !>   Compute the roots of a given Legendre polynomial. The roots are used
    !!   for the function evaluation nodes for Gauss-Legendre quadrature.
    !!
    !!   @author  Nathan A. Wukie
    !!
    !!   @param[in]  nnodes   number of Gauss-Legendre quadrature nodes requested
    !!   @param[out] nodes    Gauss-Legendre quadrature node locations
    !!
    !----------------------------------------------------------------------------------------
    subroutine gl_nodes(nnodes,nodes)
        integer(ik),    intent(in)      :: nnodes
        real(rk),       intent(inout)   :: nodes(:)

        integer(ik)                     :: i,j,polyterm
        real(rk)                        :: resid,tol,x,xnew,theta,n


        !
        ! Tolerance tested up to 10th order
        !
        !tol = 1.2e-16_rk
        tol = 1.1_rk * epsilon(1._rk)

        !
        ! compute Legendre term used to generate the correct number of nodes
        !
        polyterm = nnodes+1

        !
        ! Loop through roots so they are ordered from lowest to highest
        !
        j=1
        do i=nnodes,1,-1

            !
            ! force the condition into the loop
            !
            resid = tol + ONE


            !
            ! Use Newtons method with Chebyshev root as initial guess
            !x = cos(pi*(2._rk*real(i,rk) - 1._rk)/(2._rk*real(nnodes,rk)))
            !
            ! new initial guess, based on Lether
            ! J. Comput. Appl. Math. 4, 47 (1978)
            !
            n = real(nnodes,rk)
            theta = PI*(FOUR*real(i,rk) - ONE)/(FOUR*n + TWO)
            x = (ONE  + (n-ONE)/(EIGHT*n*n*n) - (ONE/(384._rk*n*n*n*n))*( &
                39._rk - 28._rk/(sin(theta)**(TWO))))*cos(theta)


            !
            ! Newton iteration to desired tolerance
            !
            do while (resid > tol)
                xnew = x - LegendreVal1D(polyterm,x)/DLegendreVal1D(polyterm,x)
                resid = abs(xnew-x)
                x = xnew
            end do
            nodes(j) = x
            j=j+1 ! storage location for correct ordering

        end do

    end subroutine gl_nodes
    !*****************************************************************************************











    !>  Compute the weights for Gauss-Legendre quadrature rules.
    !!
    !!  @author  Nathan A. Wukie
    !!
    !!  @param[in] nweights     number of Gauss-Legendre quadrature weights requested
    !!  @param[out] weights     quadrature weights associated with Gauss-Legendre nodes
    !!
    !-----------------------------------------------------------------------------------------
    subroutine gl_weights(nweights,weights)
        integer(ik),    intent(in)    :: nweights
        real(rk),       intent(inout) :: weights(:)

        real(rk)        :: nodes(nweights)
        integer(ik)     :: i,polyterm
        real(rk)        :: xi,DPoly

        !
        ! Get quadrature nodes
        !
        call gl_nodes(nweights,nodes)

        
        !
        ! compute legendre term used to generate the correct number of roots
        !
        polyterm = nweights+1



        do i=1,nweights

            !
            ! Get quadrature node and polynomial derivative
            !
            xi = nodes(i)
            DPoly = DLegendreVal1D(polyterm,xi)

            !
            ! Compute weight
            !
            weights(i) = TWO/((ONE - xi*xi)*(DPoly*DPoly))

        end do


    end subroutine gl_weights
    !******************************************************************************************












end module mod_gauss_legendre
