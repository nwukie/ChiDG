module fcn_monopole
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, ONE, TWO, THREE, FOUR, PI, EIGHT
    use type_point,     only: point_t
    use type_point_ad,  only: point_ad_t
    use type_function,  only: function_t
    use DNAD_D
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

        ! Set function name
        call self%set_name("monopole")

        ! Set function options to default settings
        call self%add_option('a',1._rk)

    end subroutine init
    !*************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!  Removed the use of complex numbers with algebraic operations to allow the use
    !!  of AD_D type for coordinates.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   8/19/2018
    !!
    !----------------------------------------------------------------------------------
    impure elemental function compute(self,time,coord) result(val)
        class(monopole_f),  intent(inout)   :: self
        real(rk),           intent(in)      :: time
        type(point_ad_t),   intent(in)      :: coord

        type(AD_D)  :: val

        type(AD_D)  :: x, y, z, r, B_real, B_imag, e_real, e_imag,  &
                       exp_real, exp_imag, mul_real, mul_imag
        real(rk)    :: omega, k, q, imag_real, imag_imag,           &
                       rhobar, cbar, a, Lw, freq

        complex(rk) :: B, imag


        !imag = cmplx(ZERO,ONE)
        imag_real = ZERO
        imag_imag = ONE

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
        ! (x + yi)*u = xu + yui 
        B_real = imag_real * ( rhobar*cbar*q*k/(FOUR*PI*r) )
        B_imag = imag_imag * ( rhobar*cbar*q*k/(FOUR*PI*r) )
        !B = imag*rhobar*cbar*q*k/(FOUR*PI*r)




        !
        ! Get indicator to compute 'real' or 'imaginary' component
        !
        a   = self%get_option_value('a')

        
        !
        ! Decompose complex number calculations into algebraic steps
        !
        ! -imag*k*r => -(x + yi)*u = -xu - yui 
        e_real = -imag_real * k * r
        e_imag = -imag_imag * k * r
        ! exp(-imag*k*r) => e^(x+yi) = e^(x) * e^(yi)
        !                => e^(yi)   = cos(y) + isin(y)
        !                => e^(x+yi) = e^(x) * cos(y) + (e^(x) * sin(y))i
        exp_real = exp(e_real) * dcos(e_imag)
        exp_imag = exp(e_real) * dsin(e_imag)
        ! B * exp(-imag*k*r) => (x + yi)*(u + vi) = (xu - yv) + (xv + yu)i
        mul_real = B_real * exp_real - B_imag * exp_imag
        mul_imag = B_real * exp_imag + B_imag * exp_real

        if ( a < ONE ) then
            !val = B*cos(-k*r)
            !val = real(B*exp(-imag*k*r))
            val = mul_real
        else if ( a > ONE ) then
            !val = B*sin(-k*r)
            !val = aimag(B*exp(-imag*k*r))
            val = mul_imag
        end if

    end function compute
    !***********************************************************************************


end module fcn_monopole
