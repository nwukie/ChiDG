module mod_blending_functions 
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE, FOUR, FIVE, EIGHT, PI
    implicit none


contains

    function transition_polynomial_order5(s) result(b)
        real(rk) :: s
        real(rk) :: b

        b = 10.0_rk*s**THREE - 15.0_rk*s**FOUR + 6.0_rk*s**FIVE

    end function transition_polynomial_order5

end module mod_blending_functions
