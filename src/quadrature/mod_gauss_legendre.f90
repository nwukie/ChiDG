!>  Gauss-Legendre quadrature rule.
!!
!!  Procedures:
!!  ------------------------------
!!  quadrature_nodes        ! general rule, needs (nterms, level, (1,2,3)-dim)
!!  quadrature_weights      ! general rule, needs (nterms, level, (1,2,3)-dim)
!!
!!  gl_nodes                ! 1D rule, only needs nnodes
!!  gl_weights              ! 1D rule, only needs nnodes
!!
!!  @author Nathan A. Wukie
!!
!---------------------------------------------------------------------------------------------
module mod_gauss_legendre
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: PI, ONE, TWO, THREE, FOUR, EIGHT, ZERO
    use mod_legendre,       only: legendre_val1D, dlegendre_val1D
    implicit none

contains








    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/28/2017
    !!
    !------------------------------------------------------------------------
    function quadrature_nodes(nterms,level,dim) result(nodes_)
        integer(ik),    intent(in)  :: nterms
        integer(ik),    intent(in)  :: level
        integer(ik),    intent(in)  :: dim

        ! Level for quadrature nodes corresponds to the number of nodes
        ! in the 1D tensor construction of the node set.
        real(rk),   allocatable, dimension(:)   :: xi_nodes, eta_nodes, zeta_nodes
        real(rk),   allocatable, dimension(:,:) :: nodes_
        integer(ik)     :: ixi, ieta, izeta, inode, nnodes, &
                           ierr, nterms1d, nnodes1d,        &
                           nnodes_xi, nnodes_eta, nnodes_zeta


        !
        ! Compute number of terms in 1D polynomial expansion
        !
        nterms1d = 0
        do while (nterms1d*nterms1d*nterms1d < nterms)
            nterms1d = nterms1d + 1
        end do


        !
        ! Compute number of nodes in quadrature set
        !
        select case (level)
            case(1)
                nnodes1d = nterms1d         ! Collocation quadrature
            case(2)
                nnodes1d = ceiling(3._rk*real(nterms1d,rk)/2._rk)
            case(3)
                nnodes1d = 2*nterms1d + 1
            case(4)
                nnodes1d = 3*nterms1d + 1
            case(5)
                nnodes1d = 5*nterms1d + 1
            case default
                call chidg_signal(FATAL, "compute_nnodes_integration: Value for gq_rule, specifying the rule for selecting number of quadrature points was not valid. Recognized values are gq_rule = (1, 2, 3)")
        end select



        !
        ! Allocate 1d arrays
        !
        allocate(xi_nodes(nnodes1d), eta_nodes(nnodes1d), zeta_nodes(nnodes1d), stat=ierr)
        if (ierr /= 0) call AllocationError




        !
        ! Get 1D nodes
        !
        select case(dim)
            case(1)
                nnodes_xi   = nnodes1d
                nnodes_eta  = 1
                nnodes_zeta = 1
                xi_nodes = gl_nodes(nnodes1d)
                eta_nodes  = ZERO
                zeta_nodes = ZERO
            case(2)
                nnodes_xi   = nnodes1d
                nnodes_eta  = nnodes1d
                nnodes_zeta = 1
                xi_nodes  = gl_nodes(nnodes1d)
                eta_nodes = gl_nodes(nnodes1d)
                zeta_nodes = ZERO
            case(3)
                nnodes_xi   = nnodes1d
                nnodes_eta  = nnodes1d
                nnodes_zeta = nnodes1d
                xi_nodes   = gl_nodes(nnodes1d)
                eta_nodes  = gl_nodes(nnodes1d)
                zeta_nodes = gl_nodes(nnodes1d)
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
    function quadrature_weights(nterms,level,dim) result(weights)
        integer(ik),    intent(in)  :: nterms
        integer(ik),    intent(in)  :: level
        integer(ik),    intent(in)  :: dim

        ! Level for quadrature weights corresponds to the number of weights
        ! in the 1D tensor construction of the node set.
        real(rk),   allocatable, dimension(:)   :: xi_weights, eta_weights, zeta_weights, weights
        integer(ik)                             :: ixi, ieta, izeta, inode, nweights,   &
                                                   ierr, nnodes1d, nterms1d,            &
                                                   nweights_xi, nweights_eta, nweights_zeta


        !
        ! Compute number of terms in 1D polynomial expansion
        !
        nterms1d = 0
        do while (nterms1d*nterms1d*nterms1d < nterms)
            nterms1d = nterms1d + 1
        end do


        !
        ! Compute number of nodes in quadrature set
        !
        select case (level)
            case(1)
                nnodes1d = nterms1d         ! Collocation quadrature
            case(2)
                nnodes1d = ceiling(3._rk*real(nterms1d,rk)/2._rk)
            case(3)
                nnodes1d = 2*nterms1d + 1
            case(4)
                nnodes1d = 3*nterms1d + 1
            case(5)
                nnodes1d = 5*nterms1d + 1
            case default
                call chidg_signal(FATAL, "compute_nnodes_integration: Value for gq_rule, specifying the rule for selecting number of quadrature points was not valid. Recognized values are gq_rule = (1, 2, 3)")
        end select







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
                xi_weights = gl_weights(nnodes1d)
                eta_weights  = ONE
                zeta_weights = ONE
            case(2)
                nweights_xi   = nnodes1d
                nweights_eta  = nnodes1d
                nweights_zeta = 1
                xi_weights   = gl_weights(nnodes1d)
                eta_weights  = gl_weights(nnodes1d)
                zeta_weights = ONE
            case(3)
                nweights_xi   = nnodes1d
                nweights_eta  = nnodes1d
                nweights_zeta = nnodes1d
                xi_weights   = gl_weights(nnodes1d)
                eta_weights  = gl_weights(nnodes1d)
                zeta_weights = gl_weights(nnodes1d)
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
    function gl_nodes(nnodes) result(nodes)
        integer(ik),    intent(in)      :: nnodes

        real(rk)    :: nodes(nnodes)
        integer(ik) :: i,j,polyterm
        real(rk)    :: resid,tol,x,xnew,theta,n


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
                xnew = x - legendre_val1D(polyterm,x)/dlegendre_val1D(polyterm,x)
                resid = abs(xnew-x)
                x = xnew
            end do
            nodes(j) = x
            j=j+1 ! storage location for correct ordering

        end do

    end function gl_nodes
    !*****************************************************************************************











    !>  Compute the weights for Gauss-Legendre quadrature rules.
    !!
    !!  @author  Nathan A. Wukie
    !!
    !!  @param[in] nweights     number of Gauss-Legendre quadrature weights requested
    !!  @param[out] weights     quadrature weights associated with Gauss-Legendre nodes
    !!
    !-----------------------------------------------------------------------------------------
    function gl_weights(nweights) result(weights)
        integer(ik),    intent(in)    :: nweights

        real(rk)        :: nodes(nweights), weights(nweights)
        integer(ik)     :: i,polyterm
        real(rk)        :: xi,DPoly

        !
        ! Get quadrature nodes
        !
        nodes = gl_nodes(nweights)

        
        !
        ! compute legendre term used to generate the correct number of roots
        !
        polyterm = nweights+1



        do i=1,nweights

            ! Get quadrature node and polynomial derivative
            xi = nodes(i)
            DPoly = dlegendre_val1D(polyterm,xi)

            ! Compute weight
            weights(i) = TWO/((ONE - xi*xi)*(DPoly*DPoly))

        end do


    end function gl_weights
    !******************************************************************************************












end module mod_gauss_legendre
