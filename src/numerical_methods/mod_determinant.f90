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








end module mod_determinant
