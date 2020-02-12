module fcn_scalar_adv_diff_bl_solution
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO, ONE, TWO, THREE, FIVE
    use type_point,     only: point_t
    use type_point_ad,  only: point_ad_t
    use type_function,  only: function_t
    use DNAD_D
    implicit none





    !>  scalar_adv_diff_bl_solution function.
    !!
    !!  Three 1D scalar_adv_diff_bl_solution functions are computed; one for each dimension. They are multiplied
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
    type, extends(function_t), public :: scalar_adv_diff_bl_solution_f


    contains

        procedure   :: init
        procedure   :: compute

    end type scalar_adv_diff_bl_solution_f
    !*************************************************************************************



contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !-------------------------------------------------------------------------
    subroutine init(self)
        class(scalar_adv_diff_bl_solution_f),  intent(inout)   :: self

        !
        ! Set function name
        !
        call self%set_name("scalar_adv_diff_bl_solution")


        !
        ! Set function options to default settings
        !
        call self%add_option('mu',0.05_rk)
        call self%add_option('cx',1._rk)


    end subroutine init
    !*************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !----------------------------------------------------------------------------------
    impure elemental function compute(self,time,coord) result(val)
        class(scalar_adv_diff_bl_solution_f),  intent(inout)  :: self
        real(rk),           intent(in)  :: time
        type(point_ad_t),   intent(in)  :: coord

        type(AD_D)                      :: val

        
        integer(ik) :: fcn_dim
        type(AD_D)  :: x, y, z
        real(rk)    :: mu, cx 

        !
        ! Get inputs and function parameters
        !
        x = coord%c1_
        y = coord%c2_
        z = coord%c3_

        mu = self%get_option_value('mu')
        cx = self%get_option_value('cx')



        val = (ONE-exp((x-ONE)*(cx/mu)))/(ONE-exp(-cx/mu))
    end function compute
    !***********************************************************************************


end module fcn_scalar_adv_diff_bl_solution
