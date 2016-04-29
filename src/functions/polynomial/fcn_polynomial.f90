module fcn_polynomial
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX
    use type_function,  only: function_t
    use type_point,     only: point_t
    implicit none
    private



    !>  A polynomial function: up to sixth-order
    !!
    !!  \f$   f(x) = ax^6 + bx^5 + cx^4 + dx^3 + ex^2 + fx + g  \f$
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/14/2016
    !!
    !-------------------------------------------------------------------------
    type, extends(function_t), public :: polynomial_f

    contains

        procedure   :: init
        procedure   :: compute

    end type polynomial_f
    !**************************************************************************



contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(polynomial_f),  intent(inout)   :: self

        !
        ! Set function name
        !
        call self%add_name("polynomial")


        !
        ! Set function options to default settings
        !
        ! Options here are polynomial coefficients:
        !   f(x) = ax^6 + bx^5 + cx^4 + dx^3 + ex^2 + fx + g
        !
        call self%add_option('a', ZERO)
        call self%add_option('b', ZERO)
        call self%add_option('c', ZERO)
        call self%add_option('d', ZERO)
        call self%add_option('e', ZERO)
        call self%add_option('f', ZERO)
        call self%add_option('g', ZERO)


    end subroutine init
    !*************************************************************************




    !>
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!
    !--------------------------------------------------------------------------
    impure elemental function compute(self,time,coord) result(val)
        class(polynomial_f),    intent(inout)   :: self
        real(rk),               intent(in)      :: time
        type(point_t),          intent(in)      :: coord

        real(rk)    :: a, b, c, d, e, f, g
        real(rk)    :: x
        real(rk)    :: val

        ! Get x
        x = coord%c1_

        ! Get function options
        a = self%get_option_value('a')
        b = self%get_option_value('b')
        c = self%get_option_value('c')
        d = self%get_option_value('d')
        e = self%get_option_value('e')
        f = self%get_option_value('f')
        g = self%get_option_value('g')


        !   f(x) = ax^6 + bx^5 + cx^4 + dx^3 + ex^2 + fx + g
        val = a*x**SIX + b*x**FIVE + c*x**FOUR + d*x**THREE + e*x**TWO + f*x + g

    end function compute
    !**************************************************************************


end module fcn_polynomial
