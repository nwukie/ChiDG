module fcn_ysquared
    use mod_kinds,      only: rk,ik
    use type_function,  only: function_t
    use type_point,     only: point_t
    implicit none
    private





    !> y-squared function.
    !!
    !!  \f$     f(t,\vec{x}) = y^2      \f$
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------------------
    type, extends(function_t), public :: ysquared_f

    contains

        procedure   :: order
        procedure   :: calc
        
    end type ysquared_f
    !************************************************************************************



contains


    function order(self)
        class(ysquared_f), intent(in)  :: self
        integer(ik)                 :: order

        order = 3

    end function



    elemental function calc(self,pt)
        class(ysquared_f),  intent(in)  :: self
        type(point_t),      intent(in)  :: pt
        real(rk)                        :: calc

        ! f(x) = y**2
        calc = pt%c2_  *  pt%c2_

    end function


end module fcn_ysquared
