module fcn_gaussian
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, ONE, TWO, THREE, FIVE
    use type_point,     only: point_t
    use atype_function, only: function_t
    implicit none
    private

    type, extends(function_t), public :: gaussian_f
        private


        ! constants in the gaussian function
        !
        ! f(x) = a exp(- (x-b)**2 / 2c**2)
        !

        real(rk)    :: a = ONE
        real(rk)    :: b_x = ZERO
        real(rk)    :: b_y = ZERO
        real(rk)    :: b_z = ZERO

        !real(rk)    :: c = THREE
        real(rk)    :: c = ONE

    contains
        procedure   :: order
        procedure   :: calc
    end type gaussian_f



contains


    function order(self)
        class(gaussian_f), intent(in)   :: self
        integer(ik)                     :: order

        order = 3

    end function



    elemental function calc(self,pt)
        class(gaussian_f),  intent(in)  :: self
        type(point_t),      intent(in)  :: pt
        real(rk)                        :: calc

        real(rk)    :: x,   y,   z, &
                       v_x, v_y, v_z

        x = pt%c1_
        y = pt%c2_
        z = pt%c3_

        v_x = self%a * exp( - ((x - self%b_x)**TWO) / (TWO * self%c**TWO))
        v_y = self%a * exp( - ((y - self%b_y)**TWO) / (TWO * self%c**TWO))
        v_z = self%a * exp( - ((z - self%b_z)**TWO) / (TWO * self%c**TWO))


        calc = v_x * v_y * v_z

    end function


end module fcn_gaussian
