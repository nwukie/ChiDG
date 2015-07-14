module mod_testfcn
    use mod_kinds,      only: rk,ik
    use atype_function, only: function_t
    implicit none
    private

    type, extends(function_t), public :: testfcn

    contains
        procedure   :: order
        procedure   :: calc
    end type testfcn



contains


    function order(self)
        class(testfcn), intent(in)  :: self
        integer(ik)                 :: order

        order = 3

    end function



    elemental function calc(self,pt)
        class(testfcn), intent(in)  :: self
        type(point_t),  intent(in)  :: pt
        real(rk)                    :: calc

        calc = 1.0_rk

    end function


end module mod_testfcn
