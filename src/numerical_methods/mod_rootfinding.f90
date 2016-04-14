module mod_rootfinding
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: ZERO, ONE, TWO
    use type_point,     only: point_t
    use type_function,  only: function_t
    implicit none





contains




    !>  Use the bisection method to find the root of a function, given
    !!  two locations that are known to bound the root. 
    !!
    !!  Note: because functions are defined in ChiDG as functions of space and time, we
    !!        set time here to zero, and use the first coordinate of the point data type
    !!        point%c1_ for our working coordinate.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/13/2016
    !!
    !!  @param[in]  fcn         Function definition that can be evaluated.
    !!  @param[in]  xlower      Lower location, smaller than the location of the root.
    !!  @param[in]  xupper      Upper location, larger than the location of the root.
    !!  @param[in]  tol         Optional tolerance on convergnance
    !!
    !------------------------------------------------------------------------------------
    function bisect(fcn, xlower, xupper, tol) result(fc)
        class(function_t),  intent(inout)           :: fcn
        real(rk),           intent(in)              :: xlower
        real(rk),           intent(in)              :: xupper
        real(rk),           intent(inout), optional :: tol

        logical         :: not_zero = .true.
        real(rk)        :: fa, fb, fc, time
        type(point_t)   :: a, b, c

        call a%set(ZERO, ZERO, ZERO)
        call b%set(ZERO, ZERO, ZERO)
        call c%set(ZERO, ZERO, ZERO)
        time = ZERO


        !
        ! Set default tolerance
        !
        if ( .not. present(tol) ) then
            tol = 1.e-14_rk
        end if


        !
        ! Initialize bounds
        !
        a%c1_ = xlower    ! lower bound - 'a'
        b%c1_ = xupper    ! upper bound - 'b'
        fa = fcn%compute(time,a) ! function value at 'a'
        fb = fcn%compute(time,b) ! function value at 'b'


        !
        ! Bisection iteration
        !
        do while ( not_zero ) 

            ! Compute new location
            c%c1_ = (a%c1_ + b%c1_)/TWO

            ! Compute new function value at location 'c'
            fc = fcn%compute(time,c)

            ! Check for convergence
            not_zero = ( fc > tol ) .or. ((b%c1_-a%c1_)/TWO > tol)

            ! Update bounds
            if ( sign(fc,fc) == sign(fa,fa) ) then
                a = c
                fa = fcn%compute(time,a)
            else
                b = c
                fb = fcn%compute(time,b)
            end if

        end do


    end function bisect
    !*************************************************************************************











end module mod_rootfinding
