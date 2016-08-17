module mod_condition
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ONE

    use mod_inv,    only: inv

    implicit none

    ! External procedures defined in LAPACK
    external DGETRI
    external DGETRF
    external DGECON
    external DLANGE
    real(rk) :: DLANGE

contains

!    subroutine condition_number(dim1,dim2,mat,rcond)
!        integer(kind=ik), intent(in)                        :: dim1,dim2
!        real(kind=rk),    intent(in), dimension(dim1,dim2)  :: mat
!        real(kind=rk),    intent(inout)                     :: rcond
!
!
!        real(kind=rk)                                   :: mat_norm,invmat_norm
!        integer(kind=ik)                                :: i, info, j, n, idim1,idim2
!        real(kind=rk),      dimension(dim1,dim2)        :: invmat
!        real(kind=rk),      dimension(dim2)             :: vector_sum_mat, vector_sum_invmat
!
!        invmat = inv(mat)
!
!        ! compute the 1-norm of each matrix
!        do idim2 = 1,dim2
!            vector_sum_mat(idim2)    = sum(abs(mat(:,idim2)))
!            vector_sum_invmat(idim2) = sum(abs(invmat(:,idim2)))
!        end do
!
!        mat_norm    = maxval(vector_sum_mat)
!        invmat_norm = maxval(vector_sum_invmat)
!
!        ! compute the matrix condition number
!        rcond = mat_norm*invmat_norm
!
!    end subroutine condition_number
!





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!  TODO: Try to find URL for citation.
    !!
    !---------------------------------------------------------------------------
    function cond(A) result(rcond_out)
        real(rk), dimension(:,:), intent(in) :: A
        real(rk), dimension(size(A,1),size(A,2)) :: Ainv

!        real(rk), dimension(size(A,1)) :: work    ! work array for LAPACK
        real(rk), dimension(size(A,1)) :: nwork
        real(rk), dimension(size(A,1)*4) :: conwork
        integer,  dimension(size(A,1))   :: iconwork
        integer,  dimension(size(A,1)) :: ipiv    ! pivot indices
        integer :: n, info

        real(rk)    :: anorm, rcond
        real(rk)    :: rcond_out



        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)


        !
        ! Compute matrix 1-norm
        !
        anorm = DLANGE( 'I', n, n, A, n, nwork)
      

        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call DGETRF(n, n, Ainv, n, ipiv, info)

        if (info /= 0) then
            stop 'Matrix is numerically singular!'
        end if



        !
        ! DGECON Estimate condition number
        !
        call DGECON('I',n,Ainv,n,anorm,RCOND,conwork,iconwork, info)
        !print*, 'Jacobi condition number: ', ONE/RCOND

        !
        ! RCOND is reciprical of condition number, so invert for output
        !
        rcond_out = ONE/rcond


    end function cond
    !******************************************************************************


















end module mod_condition
