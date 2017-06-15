!>  Legendre-Gauss quadrature rule.
!!
!!  Procedures:
!!  ------------------------------
!!  lg_nodes
!!  lg_weights
!!
!!  @author Nathan A. Wukie
!!
!---------------------------------------------------------------------------------------------
module mod_legendre_gauss
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: PI, ONE, TWO, THREE, FOUR, EIGHT
    use mod_legendre,       only: LegendreVal1D, DLegendreVal1D
    implicit none

contains





    !>   Compute the roots of a given Legendre polynomial. The roots are used
    !!   for the function evaluation nodes for Legendre-Gauss quadrature.
    !!
    !!   @author  Nathan A. Wukie
    !!
    !!   @param[in]  res    number of Legendre-Gauss quadrature nodes requested.
    !!   @param[out] nodes  Legendre-Gauss quadrature node locations.
    !!
    !----------------------------------------------------------------------------------------
    function lg_nodes(res) result(nodes_)
        integer(ik),    intent(in)  :: res

        real(rk)    :: nodes_(res)

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
        polyterm = res+1

        !
        ! Loop through roots so they are ordered from lowest to highest
        !
        j=1
        do i=res,1,-1

            !
            ! force the condition into the loop
            !
            resid = tol + ONE


            !
            ! Use Newtons method with Chebyshev root as initial guess
            !x = cos(pi*(2._rk*real(i,rk) - 1._rk)/(2._rk*real(res,rk)))
            !
            ! new initial guess, based on Lether
            ! J. Comput. Appl. Math. 4, 47 (1978)
            !
            n = real(res,rk)
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
            nodes_(j) = x
            j=j+1 ! storage location for correct ordering

        end do

    end function lg_nodes
    !*****************************************************************************************











    !>  Compute the weights for Gauss-Legendre quadrature rules.
    !!
    !!  @author  Nathan A. Wukie
    !!
    !!  @param[in]  res         number of Legendre-Gauss quadrature weights requested
    !!  @param[out] weights_    quadrature weights associated with Legendre-Gauss nodes
    !!
    !-----------------------------------------------------------------------------------------
    function lg_weights(res) result(weights_)
        integer(ik),    intent(in)    :: res

        real(rk)    :: weights_(res), nodes_(res)
        integer(ik) :: i,polyterm
        real(rk)    :: xi,DPoly

        !
        ! Get quadrature nodes
        !
        nodes_ = lg_nodes(res)

        
        !
        ! compute legendre term used to generate the correct number of roots
        !
        polyterm = res+1



        do i=1,res

            !
            ! Get quadrature node and polynomial derivative
            !
            xi = nodes_(i)
            DPoly = DLegendreVal1D(polyterm,xi)

            !
            ! Compute weight
            !
            weights_(i) = TWO/((ONE - xi*xi)*(DPoly*DPoly))

        end do


    end function lg_weights
    !******************************************************************************************












end module mod_legendre_gauss
