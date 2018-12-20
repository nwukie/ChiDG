module model_mnp_shock_sensor
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, TWO, PI
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    use ieee_arithmetic 
    implicit none


    


    !> Int. J. Numer. Meth. Fluids 2016; 82:398–416
    !! Dilation-based shock capturing for high-order methods
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    09/07/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(model_t)  :: mnp_shock_sensor_t


    contains

        procedure   :: init
        procedure   :: compute

    end type mnp_shock_sensor_t
    !***************************************************************************************





contains




    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    09/07/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)   
        class(mnp_shock_sensor_t), intent(inout)   :: self

        call self%set_name('MNP Shock Sensor')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('MNP Shock Sensor')
    end subroutine init
    !***************************************************************************************






    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    09/07/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(mnp_shock_sensor_t),  intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        integer(ik) :: order
        type(AD_D), allocatable,    dimension(:) ::         &
            div_vel, D, D_smooth, c_star, temp1, temp2

        real(rk),   allocatable,    dimension(:) :: h_ref
        real(rk) :: alpha, beta

        div_vel = worker%get_field('Velocity Divergence', 'value')
        c_star = worker%get_field('Critical Sound Speed', 'value')
        h_ref = worker%element_size('interior')
        order = worker%solution_order('interior')

        if (order == 0 ) order = 1

        D = -(1.5_rk*minval(h_ref)/real(order, rk))*div_vel/c_star
        if (any(.not. (ieee_is_finite(D(:)%x_ad_)))) print *, 'D is infinite'
        if (any(.not. (ieee_is_finite(D(:)%x_ad_)))) print *, 'c_star', c_star(:)%x_ad_

        ! The reference gives the value alpha = 1.0e4, but this has resulted
        ! in values of Infinity for the exponential below.
        ! Use a smaller value of alpha, or a different soft-max variant?
        alpha = 1.0e2_rk !1.0e4_rk
        beta = 0.01_rk

        D_smooth = D
        D_smooth = log(ONE +  exp(alpha*(D-beta)))/alpha
        if (any(.not. (ieee_is_finite(D_smooth(:)%x_ad_)))) D_smooth = max(beta, D)
        !if (any(.not. (ieee_is_finite(D_smooth(:)%x_ad_)))) print *, 'D_smooth is infinite'
        !if (any(.not. (ieee_is_finite(D_smooth(:)%x_ad_)))) print *, 'D', D(:)%x_ad_
        !print *, 'D_smooth', D_smooth(:)%x_ad_

        call worker%store_model_field('MNP Shock Sensor', 'value', D_smooth)

    end subroutine compute
    !***************************************************************************************




end module model_mnp_shock_sensor
