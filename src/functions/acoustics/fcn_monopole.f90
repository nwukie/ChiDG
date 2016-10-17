module fcn_monopole
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, ONE, TWO, THREE, FOUR, PI, EIGHT
    use type_point,     only: point_t
    use type_function,  only: function_t
    implicit none





    !>  Gaussian function.
    !!
    !!  Three 1D Gaussian functions are computed; one for each dimension. They are multiplied
    !!  together to create a 3D version of the function.
    !!
    !!  \f$  f_x(t,\vec{x}) = a e^{- \frac{(x-b_x)^2}{2c^2} }    \\
    !!       f_y(t,\vec{x}) = a e^{- \frac{(y-b_y)^2}{2c^2} }    \\
    !!       f_z(t,\vec{x}) = a e^{- \frac{(z-b_z)^2}{2c^2} }    \\
    !!       f(t,\vec{x}) = f_x * f_y * f_z                      \f$
    !!
    !!  Function parameters:
    !!
    !!  \f$ a   -   \text{Amplitude of the distribution}   \\
    !!      b_i -   \text{Offset in coordinate 'i'}        \\
    !!      c   -   \text{Width of the distribution}   \f$
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-------------------------------------------------------------------------------------
    type, extends(function_t), public :: monopole_f


    contains

        procedure   :: init
        procedure   :: compute

    end type monopole_f
    !*************************************************************************************



contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(monopole_f),  intent(inout)   :: self

        !
        ! Set function name
        !
        call self%set_name("monopole")


        !
        ! Set function options to default settings
        !
        call self%add_option('a',1._rk)


    end subroutine init
    !*************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !----------------------------------------------------------------------------------
    impure elemental function compute(self,time,coord) result(val)
        class(monopole_f),  intent(inout)  :: self
        real(rk),           intent(in)  :: time
        type(point_t),      intent(in)  :: coord

        real(rk)                        :: val

        real(rk)    :: x,   y,   z,   r,    &
                       omega, k, q,         &
                       rhobar, cbar, a, Lw, freq

        complex(rk) :: B, imag


        imag = cmplx(ZERO,ONE)

        !
        ! Get inputs and function parameters
        !
        x = coord%c1_
        y = coord%c2_
        z = coord%c3_
        r = sqrt(x*x + y*y + z*z)


        !
        ! Set physical constants
        !
        rhobar = 1.225_rk
        cbar   = 340.25_rk
        Lw     = 0.004819_rk
        freq   = 42.5_rk


        !
        ! Compute monpole amplitude(B) at r
        !
        omega = TWO*PI*freq
        k     = omega/cbar
        q = sqrt(EIGHT*PI*cbar*Lw/(rhobar*((TWO*PI*freq)**TWO)))
        B = imag*rhobar*cbar*q*k/(FOUR*PI*r)




        !
        ! Get indicator to compute 'real' or 'imaginary' component
        !
        a   = self%get_option_value('a')


        if ( a < ONE ) then
            !val = B*cos(-k*r)
            val = real(B*exp(-imag*k*r))
        else if ( a > ONE ) then
            !val = B*sin(-k*r)
            val = aimag(B*exp(-imag*k*r))
        end if

    end function compute
    !***********************************************************************************


end module fcn_monopole
