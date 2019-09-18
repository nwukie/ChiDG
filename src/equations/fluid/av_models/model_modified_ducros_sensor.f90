module model_modified_ducros_sensor
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, TWO, PI
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
    type, extends(model_t)  :: modified_ducros_sensor_t


    contains

        procedure   :: init
        procedure   :: compute

    end type modified_ducros_sensor_t
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
        class(modified_ducros_sensor_t), intent(inout)   :: self

        call self%set_name('Modified Ducros Sensor')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Shock Sensor')
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
        class(modified_ducros_sensor_t),  intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        integer(ik) :: ii
        type(AD_D), allocatable,    dimension(:) ::         &
            grad1_u,       grad2_u,       grad3_u,          &
            grad1_v,       grad2_v,       grad3_v,          &
            grad1_w,       grad2_w,       grad3_w,          &
            curl_vel_1, curl_vel_2, curl_vel_3,             &
            div_vel, density, mom1, mom2, mom3, vel_mag, cutoff, D, Dm

        real(rk),   allocatable,    dimension(:) :: h
        real(rk) :: Dmax, d0, delta_d

        density = worker%get_field('Density',    'value')
        mom1    = worker%get_field('Momentum-1', 'value')
        mom2    = worker%get_field('Momentum-2', 'value')
        mom3    = worker%get_field('Momentum-3', 'value')


        vel_mag = (mom1**TWO+mom2**TWO+mom3**TWO)/(density**TWO ) 
        if (any(ieee_is_nan(vel_mag(:)%x_ad_))) print*, 'vel_mag squared is nan'
        vel_mag = sqrt(vel_mag +  1.0e-12_rk)

        h = worker%element_size('interior')

        cutoff = 0.1_rk*vel_mag/minval(h) 

        div_vel = worker%get_field('Velocity Divergence', 'value')
        if (any(ieee_is_nan(div_vel(:)%x_ad_))) print*, 'div_vel is nan'

        curl_vel_1 = worker%get_field('Velocity Curl - 1', 'value')
        curl_vel_2 = worker%get_field('Velocity Curl - 2', 'value')
        curl_vel_3 = worker%get_field('Velocity Curl - 3', 'value')

        if (any(ieee_is_nan(curl_vel_1(:)%x_ad_))) print*, 'curl_vel1 is nan'
        if (any(ieee_is_nan(curl_vel_2(:)%x_ad_))) print*, 'curl_vel2 is nan'
        if (any(ieee_is_nan(curl_vel_3(:)%x_ad_))) print*, 'curl_vel3 is nan'
        D = div_vel**TWO/(div_vel**TWO+curl_vel_1**TWO+curl_vel_2**TWO+curl_vel_3**TWO + 1.0e-9_rk)

        Dm = div_vel**TWO/(div_vel**TWO+ cutoff**TWO + 1.0e-9_rk)
        


        !D = D*Dm
        
        !Dmax = maxval(D(:)%x_ad_)

        !D = Dmax
        do ii = 1, size(D)
            !if ((ieee_is_nan(cutoff(ii)%x_ad_))) then
            !    print *, "cutoff is nan"
            !    
            !    print *, 'cutoff', cutoff(ii)%x_ad_
            !    print *, 'vel mag', vel_mag(ii)%x_ad_
            !    print *, 'minval h', minval(h)
            !end if
            !if (any(ieee_is_nan(div_vel(ii)%xp_ad_(:)))) print *, "div vel deriv is nan"
            !if (any(ieee_is_nan(curl_vel_1(ii)%xp_ad_(:)))) print *, "curl vel 1 deriv is nan"
            !if (any(ieee_is_nan(curl_vel_2(ii)%xp_ad_(:)))) print *, "curl vel 2 deriv is nan"
            !if (any(ieee_is_nan(curl_vel_3(ii)%xp_ad_(:)))) print *, "curl vel 3 deriv is nan"
            !if (any(ieee_is_nan(D(ii)%xp_ad_(:)))) print *, "D deriv is nan"
            !if ((ieee_is_nan(D(ii)%x_ad_))) print *, "D  is nan"
            !D(ii)%xp_ad_(:) = 0.0_rk

            !if (div_vel(ii)%x_ad_ > 0.0_rk) D(ii) = 0.0_rk

            d0 = 0.5_rk
            delta_d = 0.49_rk
            if (D(ii) <= (d0-delta_d)) then

                D(ii) = 0.0_rk

            else if ((d0+delta_d) <= D(ii)) then
                    
                D(ii) = 1.0_rk
            else


                D(ii) = 0.5_rk*(1.0_rk+sin(PI*(D(ii)-d0)/(TWO*delta_d)))
            end if
        end do
        call worker%store_model_field('Shock Sensor', 'value', D)

    end subroutine compute
    !***************************************************************************************




end module model_modified_ducros_sensor
