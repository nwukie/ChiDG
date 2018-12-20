module model_fnp_shock_sensor
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, TWO, PI, RKTOL
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    use ieee_arithmetic 
    implicit none


    


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    07/11/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(model_t)  :: fnp_shock_sensor_t


    contains

        procedure   :: init
        procedure   :: compute

    end type fnp_shock_sensor_t
    !***************************************************************************************





contains




    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    07/11/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)   
        class(fnp_shock_sensor_t), intent(inout)   :: self

        call self%set_name('FNP Shock Sensor')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Shock Sensor - Beta')
        call self%add_model_field('Shock Sensor - Kappa')
        call self%add_model_field('Shock Sensor - Mu')
    end subroutine init
    !***************************************************************************************






    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    07/11/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(fnp_shock_sensor_t),  intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        integer(ik) :: ii
        type(AD_D), allocatable,    dimension(:) ::         &
            grad1_u,       grad2_u,       grad3_u,          &
            grad1_v,       grad2_v,       grad3_v,          &
            grad1_w,       grad2_w,       grad3_w,          &
            curl_vel_1, curl_vel_2, curl_vel_3,             &
            s_omega, s_d,                                   &
            s_beta, s_beta_smoothed,                        &
            s_kappa, s_kappa_smoothed,                        &
            s_mu, s_mu_smoothed,                        &
            s_thr, s_min, s_max                        &
            t_xi, t_stag, c_star, Lv_norm, v_max,   &
            div_vel, density, mom1, mom2, mom3, vel_mag, cutoff, D, Dm

        real(rk),   allocatable,    dimension(:) :: h
        real(rk) :: Dmax, d0, delta_d

        density = worker%get_field('Density',    'value')
        mom1    = worker%get_field('Momentum-1', 'value')
        mom2    = worker%get_field('Momentum-2', 'value')
        mom3    = worker%get_field('Momentum-3', 'value')


        vel_mag = (mom1**TWO+mom2**TWO+mom3**TWO)/(density**TWO ) 
        if (any(ieee_is_nan(vel_mag(:)%x_ad_))) print*, 'vel_mag squared is nan'
        vel_mag = sqrt(vel_mag +  RKTOL**TWO)

        h_ref = worker%element_size('interior')


        div_vel = worker%get_field('Velocity Divergence', 'value')
        if (any(ieee_is_nan(div_vel(:)%x_ad_))) print*, 'div_vel is nan'

        curl_vel_1 = worker%get_field('Velocity Curl - 1', 'value')
        curl_vel_2 = worker%get_field('Velocity Curl - 2', 'value')
        curl_vel_3 = worker%get_field('Velocity Curl - 3', 'value')

        if (any(ieee_is_nan(curl_vel_1(:)%x_ad_))) print*, 'curl_vel1 is nan'
        if (any(ieee_is_nan(curl_vel_2(:)%x_ad_))) print*, 'curl_vel2 is nan'
        if (any(ieee_is_nan(curl_vel_3(:)%x_ad_))) print*, 'curl_vel3 is nan'

        s_d     = div_vel**TWO/(div_vel**TWO+curl_vel_1**TWO+curl_vel_2**TWO+curl_vel_3**TWO + RKTOL**TWO)

        c_star = worker%get_field('Sound Speed at Critical Temperature', 'value')
        s_omega =  -(h_beta)/real(order, rk)*div_vel/c_star
        
        s_beta  = s_d*s_omega

        s_thr = 0.01_rk
        s_min = ZERO
        s_max = TWO/sqrt(gam**TWO-ONE)

        s_beta_smoothed = fnp_f(s_beta, s_thr, s_min, s_max)

        call worker%store_model_field('Shock Sensor - Beta', 'value', s_beta_smoothed)

        t_xi = worker%get_field('Temperature Gradient Under Reference Metric - Magnitude', 'value')
        t_stag = worker%get_field('Stagnation Temperature', 'value')

        s_kappa = (h_ref/real(order, rk))*t_xi/t_stag

        s_thr = ONE
        s_min = ZERO
        s_max = TWO

        s_kappa_smoothed = fnp_f(s_kappa, s_thr, s_min, s_max)

        call worker%store_model_field('Shock Sensor - Kappa', 'value', s_kappa_smoothed)

        Lv_norm = worker%get_field('Spectral Norm Deleted Diagonal Velocity Gradient', 'value')
        v_max   = worker%get_field('Maximum Isentropic Velocity', 'value')
        s_mu = (h_ref/real(order, rk))*Lv_norm/v_max

        s_thr = ONE
        s_min = ZERO
        s_max = TWO

        s_mu_smoothed = fnp_f(s_mu, s_thr, s_min, s_max)

        call worker%store_model_field('Shock Sensor - Mu', 'value', s_mu_smoothed)

    end subroutine compute
    !***************************************************************************************




end module model_fnp_shock_sensor
