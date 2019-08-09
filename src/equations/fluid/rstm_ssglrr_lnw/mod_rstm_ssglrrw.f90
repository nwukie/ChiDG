!> This module contains constants for the SSG-LRR-w RSTM.
!!
!! Based on the Reynolds' Stress Model from
!! https://turbmodels.larc.nasa.gov/rsm-ssglrr.html
!! Eisfeld, B., Rumsey, C., and Togiti, V., ``Verification and Validation of a Second-Moment-Closure Model," AIAA Journal, Vol. 54, No. 5, 2016, pp. 1524--1541.
!! Eisfeld, B., Rumsey, C., and Togiti, V., ``Erratum: Verification and Validation of a Second-Moment-Closure Model," AIAA Journal, Vol. 54, No. 9, 2016, p. 2926.
!!
!! @author  Eric M Wolf
!! @date    1/26/2018
!!
!!
!--------------------------------------------------------------------------------
module mod_rstm_ssglrrw
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: ONE, TWO, THREE
    implicit none



    real(rk), parameter :: SSG_LRRW_cmu             = 0.09_rk
    real(rk), parameter :: LRR_c2_parentheses       = 0.52_rk


    ! NOTE: on the TMR website, SSG = (epsilon), LRR = (omega)
    real(rk), parameter :: SSG_c1       = 1.7_rk
    real(rk), parameter :: LRR_c1       = 1.8_rk
    real(rk), parameter :: SSG_c1_star  = 0.9_rk
    real(rk), parameter :: LRR_c1_star  = 0.0_rk
    real(rk), parameter :: SSG_c2       = 1.05_rk
    real(rk), parameter :: LRR_c2       = 0.0_rk
    real(rk), parameter :: SSG_c3       = 0.8_rk
    real(rk), parameter :: LRR_c3       = 0.8_rk
    real(rk), parameter :: SSG_c3_star  = 0.65_rk
    real(rk), parameter :: LRR_c3_star  = 0.0_rk
    real(rk), parameter :: SSG_c4       = 0.625_rk
    real(rk), parameter :: LRR_c4       = 0.5_rk*(18._rk*LRR_c2_parentheses+12._rk)/11._rk
    real(rk), parameter :: SSG_c5       = 0.2_rk
    real(rk), parameter :: LRR_c5       = 0.5_rk*(-14._rk*LRR_c2_parentheses+20._rk)/11._rk


    real(rk), parameter :: SSG_D_SD       = (TWO/THREE)*0.22_rk
    real(rk), parameter :: LRR_D_SD       = 0.5_rk*SSG_LRRW_cmu


    real(rk), parameter :: SSG_D_GD       = 0.22_rk
    real(rk), parameter :: LRR_D_GD       = 0.75_rk*SSG_LRRW_cmu

    !
    ! Note: the following coefficients used to blend wall values (suffix _w) and boundary layer edge values (suffix _e)
    ! according to a blending function F1
    ! val = F1*val_w+(1-F1)*val_e
    !
    real(rk), parameter :: SSG_LRRW_alpha_e         = 0.44_rk
    real(rk), parameter :: SSG_LRRW_alpha_w         = 0.5556_rk
    real(rk), parameter :: SSG_LRRW_beta_e          = 0.0828_rk
    real(rk), parameter :: SSG_LRRW_beta_w          = 0.075_rk !Note: =0.09_rk in 1988 Wilcox k-w
    real(rk), parameter :: SSG_LRRW_sigma_e         = 0.856_rk
    real(rk), parameter :: SSG_LRRW_sigma_w         = 0.5_rk
    real(rk), parameter :: SSG_LRRW_sigma_d_e       = 1.712_rk
    real(rk), parameter :: SSG_LRRW_sigma_d_w       = 0.0_rk

<<<<<<< HEAD
    real(rk), parameter :: rstm_ssglrrw_avc        = 0.0_rk

    real(rk), parameter :: rstm_ssglrrw_k_infty        = 1.0e-3_rk
    real(rk), parameter :: rstm_ssglrrw_omega_infty        = 8538.1_rk!345.0_rk
=======
    real(rk), parameter :: rstm_ssglrrw_avc        = 75.0_rk

    real(rk), parameter :: rstm_ssglrrw_k_infty        = 4.761e-3_rk!1.0e-3_rk
    real(rk), parameter :: rstm_ssglrrw_R_infty        = (2.0_rk/3.0_rk)*rstm_ssglrrw_k_infty
    real(rk), parameter :: rstm_ssglrrw_omega_infty        = 345.0_rk!8538.1_rk!345.0_rk
>>>>>>> dev_tecio

    !
    ! Boundary condition constants
    !

    ! w_infty = rho_infty*k_infty/mu_t_ed_infty

    ! Farfield turbulence intensity
    ! k_infty = (3/2)*SSG_LRRW_tu_infty^2*U_infty^2
    real(rk), parameter :: SSG_LRRW_tu_infty    = 0.00039_rk 


    ! Farfield eddy viscosity fraction
    ! mu_t_ed_infty = SSG_LRRW_f_mu_infty*mu_infty
    real(rk), parameter :: SSG_LRRW_f_mu_infty  = 0.009_rk
end module mod_rstm_ssglrrw
