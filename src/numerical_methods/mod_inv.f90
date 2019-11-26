module mod_inv
#include <messenger.h>
    use mod_kinds,       only: rk,ik,rdouble,rsingle
    use mod_constants,   only: ONE
    use mod_determinant, only: det_3x3

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

    !>  Compute and return the inverse of a 3x3 matrix.
    !!
    !!  This implements an explicit formula for inverting a 3x3 matrix.
    !!      See: http://mathworld.wolfram.com/MatrixInverse.html 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/14/2017
    !!
    !---------------------------------------------------------------------------
    function inv_3x3(A) result(Ainv)
        real(rk),   intent(in)    :: A(3,3)

        real(rk)    :: Ainv(3,3), det

        det = det_3x3(A)

        Ainv(1,1) = ONE/det * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
        Ainv(1,2) = ONE/det * (A(1,3)*A(3,2) - A(1,2)*A(3,3))
        Ainv(1,3) = ONE/det * (A(1,2)*A(2,3) - A(1,3)*A(2,2))

        Ainv(2,1) = ONE/det * (A(2,3)*A(3,1) - A(2,1)*A(3,3))
        Ainv(2,2) = ONE/det * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
        Ainv(2,3) = ONE/det * (A(1,3)*A(2,1) - A(1,1)*A(2,3))

        Ainv(3,1) = ONE/det * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
        Ainv(3,2) = ONE/det * (A(1,2)*A(3,1) - A(1,1)*A(3,2))
        Ainv(3,3) = ONE/det * (A(1,1)*A(2,2) - A(1,2)*A(2,1))


    end function inv_3x3
    !***************************************************************************



    !>  Compute and return the derivative of the inverse of a 3x3 matrix A
    !!  given an already differentiated matrix A.
    !!
    !!  This uses a mathematical preposition for differentiating the inverse
    !!  of a square matrix
    !!      See: https://atmos.washington.edu/~dennis/MatrixCalculus.pdf 
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/17/2018
    !!
    !---------------------------------------------------------------------------
    function dinv_3x3(A,dA) result(dAinv)
        real(rk),   intent(in)    :: A(3,3)
        real(rk),   intent(in)    :: dA(3,3)

        real(rk)    :: Ainv(3,3), Ainv_n(3,3), Atemp(3,3), det, dAinv(3,3)

        ! Find determinatn of A
        !det = det_3x3(A)

        ! Find inverse of A
        Ainv = inv_3x3(A)

        ! Negate Ainv
        Ainv_n = -Ainv

        Atemp = matmul(Ainv_n,dA)

        dAinv = matmul(Atemp,Ainv)


    end function dinv_3x3
    !***************************************************************************





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



        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)


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



    !> Returns the inverse of a complex matrix calculated by finding the LU
    !! decomposition.  Depends on LAPACK.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/25/2018
    !!
    !---------------------------------------------------------------------------
    function zinv(A) result(Ainv)
        complex(kind=rk), dimension(:,:),               intent(inout)   :: A
        complex(kind=rk), dimension(size(A,1),size(A,2))                :: Ainv

        complex(kind=rk), dimension(size(A,1))   :: work    ! work array for LAPACK
        integer,          dimension(size(A,1))   :: ipiv    ! pivot indices
        integer :: n, info



        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)


        !
        ! ZGETRF computes an LU factorization of a general M-by-N complex matrix A
        ! using partial pivoting with row interchanges.
        !
        if ( rk == rdouble ) then
            call ZGETRF(n, n, Ainv, n, ipiv, info)
        else if ( rk == rsingle ) then
            call CGETRF(n, n, Ainv, n, ipiv, info)
        else
            call chidg_signal(FATAL,"inv: Invalid selected precision for matrix inversion.")
        end if

        if (info /= 0) then
            call chidg_signal(FATAL,"inv: Matrix is numerically singular!")
        end if



        !
        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        !
        if ( rk == rdouble ) then
            call ZGETRI(n, Ainv, n, ipiv, work, n, info)
        else if ( rk == rsingle ) then
            call CGETRI(n, Ainv, n, ipiv, work, n, info)
        else
            call chidg_signal(FATAL,"inv: Invalid selected precision for matrix inversion.")
        end if



        if (info /= 0) then
            call chidg_signal(FATAL,"inv: matrix inversion failed!.")
        end if

    end function zinv
    !**************************************************************************************














end module mod_inv
