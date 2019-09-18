module model_rstm_ssglrrw_lrr_coefficients
#include <messenger.h>
    use mod_kinds,              only: rk
    use mod_constants,          only: ZERO, HALF, ONE, TWO, THREE, FOUR
    use mod_rstm_ssglrrw
    use type_model,             only: model_t
    use type_chidg_worker,      only: chidg_worker_t
    use ieee_arithmetic,        only: ieee_is_nan
    use DNAD_D

    implicit none


    

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(model_t)  :: rstm_ssglrrw_lrr_coefficients_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_lrr_coefficients_t
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
        class(rstm_ssglrrw_lrr_coefficients_t), intent(inout)   :: self

        call self%set_name('RSTMSSGLRRW LRR Coefficients')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('RSTMSSGLRRW Alpha-w')
        call self%add_model_field('RSTMSSGLRRW Beta-w')
        call self%add_model_field('RSTMSSGLRRW Sigma-w')
        call self%add_model_field('RSTMSSGLRRW Sigma_d-w')

        call self%add_model_field('RSTMSSGLRRW C1')
        call self%add_model_field('RSTMSSGLRRW C1_star')
        call self%add_model_field('RSTMSSGLRRW C2')
        call self%add_model_field('RSTMSSGLRRW C3')
        call self%add_model_field('RSTMSSGLRRW C3_star')
        call self%add_model_field('RSTMSSGLRRW C4')
        call self%add_model_field('RSTMSSGLRRW C5')

        call self%add_model_field('RSTMSSGLRRW D_SD')
        call self%add_model_field('RSTMSSGLRRW D_GD')
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
        class(rstm_ssglrrw_lrr_coefficients_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        real(rk), dimension(:), allocatable :: F1, zeta
        type(AD_D), dimension(:),   allocatable ::              &
            density, distance, omega_t, omega_t_source_term, k_t, mu_l, &
            alpha_w, beta_w, sigma_w, sigma_d_w,          &
            C1, C1_star, C2, C3, C3_star, C4, C5, D_SD, D_GD

        !
        ! Interpolate solution to quadrature nodes
        !
        density             = worker%get_field('Density',                   'value') 
!        if (worker%interpolation_source == 'face exterior') then
!            distance            = worker%get_field('Wall Distance',             'value', 'face interior') 
!        else
!            distance            = worker%get_field('Wall Distance',             'value') 
!        end if
!        omega_t               = worker%get_field('Omega',                     'value')
!        omega_t_source_term   = worker%get_field('Omega Source Term',         'value')
!        k_t                 = worker%get_field('Turbulence Kinetic Energy', 'value')
!        mu_l                = worker%get_field('Laminar Viscosity',         'value')
!
!        if (any(ieee_is_nan(density(:)%x_ad_))) print *, 'density is nan'
!        if (any(ieee_is_nan(distance(:)%x_ad_))) print *, 'rstm lrr coeff distance is nan'
!        if (any(ieee_is_nan(distance(:)%x_ad_))) print *, 'rstm lrr coeff interp source ', worker%interpolation_source
!        if (any(ieee_is_nan(distance(:)%x_ad_))) then
!            if (allocated(distance)) deallocate(distance)
!            distance = density
!            distance = 1.0e-12_rk
!        end if
!        if (any(ieee_is_nan(distance(:)%x_ad_))) print *, 'distance is nan post'
!        if (any(ieee_is_nan(distance(:)%x_ad_))) print *, distance(:)%x_ad_
!        !print *, distance(:)%x_ad_
!        if (any(ieee_is_nan(omega_t(:)%x_ad_))) print *, 'omega is nan'
!        if (any(ieee_is_nan(omega_t_source_term(:)%x_ad_))) print *, 'omega source is nan'
!        if (any(ieee_is_nan(k_t(:)%x_ad_))) print *, 'k_t is nan'
!        if (any(ieee_is_nan(mu_l(:)%x_ad_))) print *, 'mu_l is nan'
!
!        zeta = density%x_ad_
!!        zeta    = min( &
!!                    max(sqrt(abs(k_t%x_ad_))/(SSG_LRRW_cmu*omega_t%x_ad_*distance%x_ad_ + 1.0e-12_rk), &
!!                    500._rk*mu_l%x_ad_/(density%x_ad_*omega_t%x_ad_*distance%x_ad_**TWO + 1.0e-12_rk)),&
!!                    FOUR*SSG_LRRW_sigma_e*density%x_ad_*k_t%x_ad_/(SSG_LRRW_sigma_d_e*omega_t_source_term%x_ad_*distance%x_ad_**TWO + 1.0e-12_rk))
!
        ! Blending coefficient
        F1 = density%x_ad_
        !F1 = tanh(zeta**FOUR)
        F1 = 1.0_rk

        if (any(ieee_is_nan(F1))) print *, 'F1 is nan'

        ! omega_t equation constants
        alpha_w     = density 
        beta_w      = density 
        sigma_w     = density 
        sigma_d_w   = density 

        alpha_w     = F1*SSG_LRRW_alpha_w   + (ONE-F1)*SSG_LRRW_alpha_e
        beta_w      = F1*SSG_LRRW_beta_w    + (ONE-F1)*SSG_LRRW_beta_e
        sigma_w     = F1*SSG_LRRW_sigma_w   + (ONE-F1)*SSG_LRRW_sigma_e
        sigma_d_w   = F1*SSG_LRRW_sigma_d_w + (ONE-F1)*SSG_LRRW_sigma_d_e

        sigma_w = 0.5_rk
        beta_w  = 0.09_rk

        call worker%store_model_field('RSTMSSGLRRW Alpha-w',    'value', alpha_w)
        call worker%store_model_field('RSTMSSGLRRW Beta-w',     'value', beta_w)
        call worker%store_model_field('RSTMSSGLRRW Sigma-w',    'value', sigma_w)
        call worker%store_model_field('RSTMSSGLRRW Sigma_d-w',  'value', sigma_d_w)

        C1      = density 
        C1_star = density 
        C2      = density 
        C3      = density 
        C3_star = density 
        C4      = density 
        C5      = density 

       ! Pressure-Strain Correlation constants
        C1      = F1*LRR_c1         + (ONE-F1)*SSG_c1
        C1_star = F1*LRR_c1_star    + (ONE-F1)*SSG_c1_star
        C2      = F1*LRR_c2         + (ONE-F1)*SSG_c2
        C3      = F1*LRR_c3         + (ONE-F1)*SSG_c3
        C3_star = F1*LRR_c3_star    + (ONE-F1)*SSG_c3_star
        C4      = F1*LRR_c4         + (ONE-F1)*SSG_c4
        C5      = F1*LRR_c5         + (ONE-F1)*SSG_c5

        call worker%store_model_field('RSTMSSGLRRW C1',         'value', C1)
        call worker%store_model_field('RSTMSSGLRRW C1_star',    'value', C1_star)
        call worker%store_model_field('RSTMSSGLRRW C2',         'value', C2)
        call worker%store_model_field('RSTMSSGLRRW C3',         'value', C3)
        call worker%store_model_field('RSTMSSGLRRW C3_star',    'value', C3_star)
        call worker%store_model_field('RSTMSSGLRRW C4',         'value', C4)
        call worker%store_model_field('RSTMSSGLRRW C5',         'value', C5)

        ! Diffusion model constants
        D_SD = density
        D_GD = density

        D_SD = F1*LRR_D_SD + (ONE-F1)*SSG_D_SD
        D_GD = F1*LRR_D_GD + (ONE-F1)*SSG_D_GD

        call worker%store_model_field('RSTMSSGLRRW D_SD', 'value', D_SD)
        call worker%store_model_field('RSTMSSGLRRW D_GD', 'value', D_GD)
    end subroutine compute
    !***************************************************************************************




end module model_rstm_ssglrrw_lrr_coefficients
