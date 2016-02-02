module fcn_xyz
    use mod_kinds,      only: rk,ik
    use type_function,  only: function_t
    use type_point,     only: point_t
    implicit none
    private





    !> xyz function
    !!
    !!  \f$     f(t,\vec{x}) = xyz      \f$
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !---------------------------------------------------------------------------------
    type, extends(function_t), public :: xyz_f

    contains

        procedure   :: order
        procedure   :: calc

    end type xyz_f
    !**********************************************************************************



contains


    function order(self)
        class(xyz_f), intent(in)  :: self
        integer(ik)               :: order

        order = 3

    end function



    elemental function calc(self,pt)
        class(xyz_f),   intent(in)  :: self
        type(point_t),  intent(in)  :: pt
        real(rk)                    :: calc

        ! f(x) = x * y * z
        calc = pt%c1_ * pt%c2_ * pt%c3_

    end function



end module fcn_xyz
