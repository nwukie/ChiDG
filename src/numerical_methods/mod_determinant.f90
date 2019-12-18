module mod_determinant
#include <messenger.h>
    use mod_kinds,      only: rk,ik,rdouble,rsingle
    use mod_constants,  only: ONE
    implicit none



contains


    !>  Compute and return the determinant of a 3x3 matrix.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/14/2017
    !!
    !!
    !----------------------------------------------------------
    function det_3x3(mat) result(det)
        real(rk),   intent(in)  :: mat(3,3)

        real(rk)    :: det

        det = mat(1,1)*mat(2,2)*mat(3,3) - mat(1,2)*mat(2,1)*mat(3,3) - &
              mat(1,1)*mat(2,3)*mat(3,2) + mat(1,3)*mat(2,1)*mat(3,2) + &
              mat(1,2)*mat(2,3)*mat(3,1) - mat(1,3)*mat(2,2)*mat(3,1)

    end function det_3x3
    !**********************************************************



    !>  Compute and return the derivative of determinant of a 3x3 matrix
    !!  wrt a given differentiated matrix
    !!
    !!  \partial det(mat) / \partial x
    !!  
    !!  dmat = dmat/dx
    !!
    !!
    !!  @author Matteo Ugolotti
    !!  @date   7/17/2018
    !!
    !!
    !----------------------------------------------------------
    function ddet_3x3(mat,dmat) result(det)
        real(rk),   intent(in)  :: mat(3,3)
        real(rk),   intent(in)  :: dmat(3,3)

        real(rk)    :: det1, det2, det
        real(rk)    :: adjugate(3,3), cofactor(3,3), temp(3,3)

        !
        ! Method 1 and 2 give the same result, 1 is less expensive
        !

        !
        ! Method 1: Direct differentiation.
        !
        !det1 =  dmat(1,1)*mat(2,2)*mat(3,3) + mat(1,1)*dmat(2,2)*mat(3,3) + mat(1,1)*mat(2,2)*dmat(3,3) &
        !      - dmat(1,2)*mat(2,1)*mat(3,3) - mat(1,2)*dmat(2,1)*mat(3,3) - mat(1,2)*mat(2,1)*dmat(3,3) &
        !      - dmat(1,1)*mat(2,3)*mat(3,2) - mat(1,1)*dmat(2,3)*mat(3,2) - mat(1,1)*mat(2,3)*dmat(3,2) &
        !      + dmat(1,3)*mat(2,1)*mat(3,2) + mat(1,3)*dmat(2,1)*mat(3,2) + mat(1,3)*mat(2,1)*dmat(3,2) &
        !      + dmat(1,2)*mat(2,3)*mat(3,1) + mat(1,2)*dmat(2,3)*mat(3,1) + mat(1,2)*mat(2,3)*dmat(3,1) &
        !      - dmat(1,3)*mat(2,2)*mat(3,1) - mat(1,3)*dmat(2,2)*mat(3,1) - mat(1,3)*mat(2,2)*dmat(3,1)
        
        !
        ! Method 2: Jacobi's formula
        !   ref: https://en.wikipedia.org/wiki/Jacobi%27s_formula
        !

        ! Find Cofactor matrix of mat
        cofactor(1,1) =  mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2)
        cofactor(1,2) = -mat(2,1)*mat(3,3)+mat(2,3)*mat(3,1)
        cofactor(1,3) =  mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1)
        cofactor(2,1) = -mat(1,2)*mat(3,3)+mat(1,3)*mat(3,2)
        cofactor(2,2) =  mat(1,1)*mat(3,3)-mat(1,3)*mat(3,1)
        cofactor(2,3) = -mat(1,1)*mat(3,2)+mat(1,2)*mat(3,1)
        cofactor(3,1) =  mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2)
        cofactor(3,2) = -mat(1,1)*mat(2,3)+mat(1,3)*mat(2,1)
        cofactor(3,3) =  mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)

        ! Adjugate of mat is the transpose of cofactor of mat
        adjugate = transpose(cofactor)

        ! Use Jacobi's formula to find the derivative of the determinant
        temp = matmul(adjugate,dmat)
        
        ! Find the trace of temp
        det2 = temp(1,1)+temp(2,2)+temp(3,3)

        ! Result
        det = det2


    end function ddet_3x3
    !**********************************************************








end module mod_determinant
