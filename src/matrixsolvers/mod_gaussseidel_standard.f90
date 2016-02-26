module mod_gaussseidel_standard
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ONE, ZERO
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


        !print*, 'A'
        !print*, A
!
!        read(*,*)






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

            
            !print*, '     ', res

            xold = x    ! store old solution




        end do ! res > self%tol




    end subroutine gaussseidel_standard



end module mod_gaussseidel_standard
