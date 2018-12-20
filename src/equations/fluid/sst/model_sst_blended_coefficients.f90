module model_sst_blended_coefficients
!#include <messenger.h>
    use mod_kinds,              only: rk
    use mod_constants,          only: ZERO, HALF, ONE, TWO, THREE, FOUR
    use type_model,             only: model_t
    use type_chidg_worker,      only: chidg_worker_t
    use ieee_arithmetic,        only: ieee_is_nan
    use mod_sst
    use DNAD_D
    implicit none


    

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(model_t)  :: sst_blended_coefficients_t

    contains

        procedure   :: init
        procedure   :: compute

    end type sst_blended_coefficients_t
    !***************************************************************************************





contains



    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)   
        class(sst_blended_coefficients_t), intent(inout)   :: self

        call self%set_name('SST Blended Coefficients')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('SST beta_k')
        call self%add_model_field('SST beta_w') 
        call self%add_model_field('SST sigma_k')
        call self%add_model_field('SST sigma_w')
        call self%add_model_field('SST sigma_d')
        call self%add_model_field('SST kappa')
        call self%add_model_field('SST alpha_w')
    

    end subroutine init
    !***************************************************************************************





    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker)
    use DNAD_D !NOTE: Somehow, without this statement, the AD sqrt function is not properly accessed, causing a compile-time error.
        class(sst_blended_coefficients_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        !real(rk), dimension(:), allocatable :: F1, zeta
        type(AD_D), dimension(:),   allocatable ::              &
            density, distance, omega_t, omega_source_term, k_t, mu_l, F1, zeta, CD, &
            beta_k, beta_w, sigma_w, sigma_k, sigma_d, kappa, alpha_w, temp1, temp2, temp3, temp4, arg

        integer(ik) :: ii

        !
        ! Interpolate solution to quadrature nodes
        !
        density             = worker%get_field('Density',                   'value') 

        if (worker%interpolation_source == 'boundary') then
            distance            = worker%get_field('Wall Distance',             'value', 'face interior') 
        else if (worker%interpolation_source == 'face exterior') then
            distance            = worker%get_field('Wall Distance',             'value', 'face interior') 
        else
            distance            = worker%get_field('Wall Distance',             'value') 
        end if

        omega_t               = worker%get_field('Omega',                     'value')
        omega_source_term   = worker%get_field('SST CD',                     'value')
        k_t                 = worker%get_field('k', 'value')

        mu_l                = worker%get_field('Laminar Viscosity',         'value')

        !print *, 'distance'
        !print *, distance(:)%x_ad_
        
        CD = density
        CD  = smax(sst_e_sigma_d*omega_source_term-1.0e-10_rk, ZERO*omega_source_term)+1.0e-10_rk 
        !CD  = max(sst_e_sigma_d*omega_source_term,1.0e-10_rk) 
        !print *, 'mu t max argument 1'
        !print *, sst_e_sigma_d*omega_source_term(:)%x_ad_


        !print *, 'CD'
        !print *, CD(:)%x_ad_
        zeta = density
        temp1 = density
        temp2 = density
        temp3 = density
        temp4 = density

        !temp1   = exp(0.5_rk*log(k_t))/(0.09_rk*exp(omega_t)*distance + 1.0e-12_rk)
        temp1   = MY_SQRT(k_t)
        !do ii = 1, size(temp1)
        !temp1(ii)%xp_ad_(:) = ZERO
        !end do
        temp1=temp1/(0.09_rk*exp(omega_t)*distance + 1.0e-14_rk)
        !print *, 'temp1'
        !print *, temp1(:)%x_ad_
        
        temp2    =  500._rk*mu_l/(density*exp(omega_t)*distance**TWO + 1.0e-14_rk)
        !print *, 'temp2'
        !print *, temp2(:)%x_ad_

        temp3    = smax(temp1, temp2)
        !print *, 'F1 -1 max argument 1'
        !print *, temp1(:)%x_ad_
        !print *, 'F1 -1 max argument 2'
        !print *, temp2(:)%x_ad_


        !temp3    = max(temp1, temp2)

        !print *, 'temp3'
        !print *, temp3(:)%x_ad_

        temp4   = TWO*sst_e_sigma_d*density*k_t/(CD*distance**TWO + 1.0e-14_rk)
        !print *, 'temp4'
        !print *, temp4(:)%x_ad_

        zeta = smin(temp3, temp4)
        !print *, 'F1-2 min argument 1'
        !print *, temp3(:)%x_ad_
        !print *, 'F1-2 min argument 2'
        !print *, temp4(:)%x_ad_


        !zeta = min(temp3, temp4)
        !print *, 'zeta'
        !print *, zeta(:)%x_ad_


        !zeta    = min( &
        !            max(sqrt((k_t(:)%x_ad_))/(0.09_rk*exp(omega_t)*distance + 1.0e-10_rk), &
        !            500._rk*mu_l/(density*exp(omega_t)*distance**TWO + 1.0e-10_rk)),&
        !            TWO*sst_e_sigma_d*density*k_t/(CD*distance**TWO + 1.0e-10_rk))

        ! Blending coefficient
        F1 = zeta
        arg = zeta**FOUR
        temp1 = F1
        temp1 = ONE-exp(-TWO*arg)
        temp2 = temp1
        temp2 = exp(-TWO*arg)+ONE
        F1 = temp1/temp2

        !do ii = 1, size(temp1)
        !    F1(ii)%xp_ad_(:) = ZERO
        !end do
        !F1 = tanh(zeta(:)%x_ad_**FOUR)
        F1 = 1.0_rk
        !F1 = 0.0_rk

        !print *, 'F1'
        !print *, F1(:)%x_ad_
        !if (any(ieee_is_nan(F1))) print *, 'F1 is nan'

        ! omega_t equation constants
        beta_k      = density
        beta_w      = density

        sigma_k     = density
        sigma_w     = density
        sigma_d     = density

        kappa       = density

        alpha_w     = density

        beta_k      = F1*sst_w_beta_k + (ONE-F1)*sst_e_beta_k
        beta_w      = F1*sst_w_beta_w + (ONE-F1)*sst_e_beta_w


        sigma_k     = F1*sst_w_sigma_k + (ONE-F1)*sst_e_sigma_k
        sigma_w     = F1*sst_w_sigma_w + (ONE-F1)*sst_e_sigma_w
        sigma_d     = F1*sst_w_sigma_d + (ONE-F1)*sst_e_sigma_d
        
        kappa       = F1*sst_w_kappa    + (ONE-F1)*sst_e_kappa

        alpha_w     = F1*sst_w_alpha_w + (ONE-F1)*sst_e_alpha_w

        call worker%store_model_field('SST beta_k',    'value', beta_k)
        call worker%store_model_field('SST beta_w',    'value', beta_w)

        call worker%store_model_field('SST sigma_k',    'value', sigma_k)
        call worker%store_model_field('SST sigma_w',    'value', sigma_w)
        call worker%store_model_field('SST sigma_d',    'value', sigma_d)

        call worker%store_model_field('SST kappa',    'value', kappa)

        call worker%store_model_field('SST alpha_w',    'value', alpha_w)
    end subroutine compute
    !***************************************************************************************




end module model_sst_blended_coefficients
