module mod_GaussLegendre
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: PI, ONE, TWO, THREE, FOUR, EIGHT
    use mod_legendre,       only: LegendreVal1D, DLegendreVal1D
    implicit none

contains





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












end module mod_GaussLegendre
