module mod_inv
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ONE

    implicit none
    ! External procedures defined in LAPACK
    external DGETRI
    external DGETRF
    external DGECON
    external DLANGE
    real(rk) :: DLANGE

contains

    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    function inv(A) result(Ainv)
      real(rk), dimension(:,:), intent(in) :: A
      real(rk), dimension(size(A,1),size(A,2)) :: Ainv

      real(rk), dimension(size(A,1)) :: work    ! work array for LAPACK
      real(rk), dimension(size(A,1)) :: nwork
      real(rk), dimension(size(A,1)*4) :: conwork
      integer,  dimension(size(A,1))   :: iconwork
      integer,  dimension(size(A,1)) :: ipiv    ! pivot indices
      integer :: n, info

      real(rk) :: anorm, rcond



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




      ! DGECON Estimate condition number
      !
      call DGECON('I',n,Ainv,n,anorm,RCOND,conwork,iconwork, info)
      !print*, 'Jacobi condition number: ', ONE/RCOND





      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call DGETRI(n, Ainv, n, ipiv, work, n, info)

      if (info /= 0) then
         stop 'Matrix inversion failed!'
      end if
    end function inv

end module mod_inv
