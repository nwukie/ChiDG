module fcn_zsquared
    use mod_kinds,      only: rk,ik
    use atype_function, only: function_t
    implicit none
    private

    type, extends(function_t), public :: zsquared_f

    contains
        procedure   :: order
        procedure   :: calc
    end type zsquared_f



contains


    function order(self)
        class(zsquared_f), intent(in)  :: self
        integer(ik)                 :: order

        order = 3

    end function



    elemental function calc(self,pt)
        class(zsquared_f),  intent(in)  :: self
        type(point_t),      intent(in)  :: pt
        real(rk)                        :: calc

        ! f(x) = z**2
        calc = pt%c3_  *  pt%c3_

    end function


end module fcn_zsquared
