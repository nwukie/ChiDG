module fcn_constant
    use mod_kinds,      only: rk,ik
    use atype_function, only: function_t
    use type_point,     only: point_t
    implicit none
    private

    type, extends(function_t), public :: constant_f

        real(rk)    :: value_   !> Constant function value

    contains
        procedure   :: order
        procedure   :: calc

        procedure   :: set
    end type constant_f



contains


    !> Set the order of the function. Should eventually be used to help
    !! determine integration rules.
    function order(self)
        class(constant_f), intent(in)  :: self
        integer(ik)                 :: order

        order = 3

    end function


    !> Function method to return function value.
    elemental function calc(self,pt)
        class(constant_f),  intent(in)  :: self
        type(point_t),      intent(in)  :: pt
        real(rk)                        :: calc

        ! f(x,y,z) = const
        calc = self%value_

    end function




    !> Function method to set constant value.
    subroutine set(self,value)
        class(constant_f),  intent(inout) :: self
        real(rk),           intent(in)    :: value

        self%value_ = value

    end subroutine


end module fcn_constant
