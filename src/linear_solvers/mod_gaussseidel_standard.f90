module mod_gaussseidel_standard
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ONE, ZERO
    use DNAD_D
!    use mod_inv,                only: inv

    implicit none
        





contains


    !> Solution routine
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!
    !!
    !--------------------------------------------------------------
    subroutine gaussseidel_standard(A,x,b)
        real(rk), intent(in)           :: A(:,:)
        real(rk), intent(inout)        :: x(:)
        real(rk), intent(in)           :: b(:)



        real(rk), allocatable, dimension(:)     :: r, diff, xold, err
        real(rk), dimension(size(x))            :: D

        integer(ik) :: iparent, ierr, iblk, ielem, i, irow, icol
        real(rk)    :: res, tol

        tol = 1.e-13_rk


        !
        ! Initialize D blocks
        !
        do irow = 1,size(A,1)
            D(irow) = ONE/A(irow,irow)
        end do


        xold = x
        x    = ZERO
        xold = ZERO



        res = 1._rk
        do while (res > tol)
            
            !
            ! form b - Ax, except for diagonal
            !
            r = b
            do irow = 1,size(A,1)
                do icol = 1,size(A,2)


                    !
                    ! Don't multiply diagonal
                    !
                    if ( irow /= icol ) then
                        r(irow) = r(irow) - A(irow,icol) * x(icol)
                    end if


                end do

                !
                ! Multiply by inverse block diagonal
                !
                x(irow) = D(irow) * r(irow)
            end do



            !
            ! Compute residual
            !
            err = b - matmul(A,x)
            res = norm2(err)
            
            xold = x    ! store old solution

        end do ! res > self%tol

    end subroutine gaussseidel_standard





    !> Solution routine
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/9/2018
    !!
    !!
    !--------------------------------------------------------------
    subroutine gaussseidel_autodiff(A,x,b)
        type(AD_D), intent(in)      :: A(:,:)
        type(AD_D), intent(inout)   :: x(:)
        type(AD_D), intent(in)      :: b(:)

        type(AD_D), allocatable, dimension(:)     :: r, diff, xold, err
        type(AD_D), dimension(size(x))            :: D
        type(AD_D)  :: res, res1

        integer(ik) :: iparent, ierr, iblk, ielem, i, irow, icol, iiter
        real(rk)    :: tol, oreduction

        tol = 1.e-11_rk


        !
        ! Initialize D blocks
        !
        do irow = 1,size(A,1)
            D(irow) = ONE/A(irow,irow)
        end do

        ! Initialize derivatives
        x = b
        x = ZERO

        xold = x
        x    = ZERO
        xold = ZERO



        res = b(1)
        res = 1._rk
        iiter = 1
        do while (res > tol)
            
            !
            ! form b - Ax, except for diagonal
            !
            r = b
            do irow = 1,size(A,1)
                do icol = 1,size(A,2)

                    !
                    ! Don't multiply diagonal
                    !
                    if ( irow /= icol ) then
                        r(irow) = r(irow) - A(irow,icol) * x(icol)
                    end if


                end do

                !
                ! Multiply by inverse block diagonal
                !
                x(irow) = D(irow) * r(irow)
            end do



            !
            ! Compute residual
            !
            err = b - matmul(A,x)
            res = norm2(err)
            if (iiter == 1) res1 = res
            
            print*, res%x_ad_
            if (log10(res/res1) < -8.) exit


            xold = x    ! store old solution


            iiter = iiter + 1
        end do ! res > self%tol

    end subroutine gaussseidel_autodiff






end module mod_gaussseidel_standard
