module mod_inv
#include <messenger.h>
    use mod_kinds,      only: rk,ik,rdouble,rsingle
    use mod_constants,  only: ONE

    implicit none
    ! External procedures defined in LAPACK
    external DGETRI
    external DGETRF
    external DGECON
    external DLANGE
    external SGETRI
    external SGETRF
    external SGECON
    external SLANGE
    real(rdouble) :: DLANGE
    real(rsingle) :: SLANGE

contains

    !> Returns the inverse of a matrix calculated by finding the LU
    !! decomposition.  Depends on LAPACK.
    !!
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------
    function inv(A) result(Ainv)
        real(rk), dimension(:,:),               intent(inout)   :: A
        real(rk), dimension(size(A,1),size(A,2))                :: Ainv

        real(rk), dimension(size(A,1))   :: work    ! work array for LAPACK
        integer,  dimension(size(A,1))   :: ipiv    ! pivot indices
        integer :: n, info

!        real(rk), dimension(size(A,1))   :: nwork
!        real(rk), dimension(size(A,1)*4) :: conwork
!        integer,  dimension(size(A,1))   :: iconwork
!        real(rk) :: anorm, rcond



        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)


!        !
!        ! Compute matrix 1-norm
!        !
!        if ( rk == rdouble ) then
!            anorm = DLANGE( 'I', n, n, A, n, nwork)
!        else if ( rk == rsingle ) then
!            anorm = SLANGE( 'I', n, n, A, n, nwork)
!        else
!            print*, "inv: error in selected precision"
!        end if



        !
        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        !
        if ( rk == rdouble ) then
            call DGETRF(n, n, Ainv, n, ipiv, info)
        else if ( rk == rsingle ) then
            call SGETRF(n, n, Ainv, n, ipiv, info)
        else
            call chidg_signal(FATAL,"inv: Invalid selected precision for matrix inversion.")
        end if

        if (info /= 0) then
            call chidg_signal(FATAL,"inv: Matrix is numerically singular!")
        end if





!        !
!        ! DGECON Estimate condition number
!        !
!        if ( rk == rdouble ) then
!            call DGECON('I',n,Ainv,n,anorm,RCOND,conwork,iconwork, info)
!        else if ( rk == rsingle ) then
!            call SGECON('I',n,Ainv,n,anorm,RCOND,conwork,iconwork, info)
!        else
!            print*, "inv: error in selected precision"
!        end if
!        !print*, 'Jacobi condition number: ', ONE/RCOND






        !
        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        !
        if ( rk == rdouble ) then
            call DGETRI(n, Ainv, n, ipiv, work, n, info)
        else if ( rk == rsingle ) then
            call SGETRI(n, Ainv, n, ipiv, work, n, info)
        else
            call chidg_signal(FATAL,"inv: Invalid selected precision for matrix inversion.")
        end if



        if (info /= 0) then
            call chidg_signal(FATAL,"inv: matrix inversion failed!.")
        end if

    end function inv
    !**************************************************************************************

end module mod_inv
