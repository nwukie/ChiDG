!> This module contains constants for the SST turbulence model as presented in
!!
!! Schoenawa, Stefan, and Ralf Hartmann. 
!! "Discontinuous Galerkin discretization of the Reynolds-averaged Navier–Stokes equations with the shear-stress transport model." 
!! Journal of Computational Physics 262 (2014): 194-216.
!!
!! @author  Eric M Wolf
!! @date    1/26/2018
!!
!!
!--------------------------------------------------------------------------------
module mod_sst
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: ONE, TWO, THREE
    implicit none

    !
    ! k-w constants
    !

    real(rk), parameter :: sst_w_beta_k             = 0.09_rk
    real(rk), parameter :: sst_w_beta_w             = 0.075_rk

    real(rk), parameter :: sst_w_sigma_k            = 0.5_rk ! BSL value
    !real(rk), parameter :: sst_w_sigma_k            = 0.85_rk ! SST value
    real(rk), parameter :: sst_w_sigma_w            = 0.5_rk
    real(rk), parameter :: sst_w_sigma_d            = 0.0_rk

    real(rk), parameter :: sst_w_kappa              = 0.41_rk

    real(rk), parameter :: sst_w_alpha_w            = sst_w_beta_w/sst_w_beta_k - sst_w_sigma_w*sst_w_kappa**TWO/sqrt(sst_w_beta_k)
    
    !
    ! k-e constants
    !

    real(rk), parameter :: sst_e_beta_k             = 0.09_rk
    real(rk), parameter :: sst_e_beta_w             = 0.0828_rk

    real(rk), parameter :: sst_e_sigma_k            = 1.0_rk
    real(rk), parameter :: sst_e_sigma_w            = 0.856_rk
    real(rk), parameter :: sst_e_sigma_d            = TWO*sst_e_sigma_w

    real(rk), parameter :: sst_e_kappa              = 0.41_rk

    real(rk), parameter :: sst_e_alpha_w              = sst_e_beta_w/sst_e_beta_k - sst_e_sigma_w*sst_e_kappa**TWO/sqrt(sst_e_beta_k)

    ! Artificial viscosity constant
    real(rk), parameter :: sst_avc                  = 100.0_rk


    !
    ! Boundary condition constants
    !

    ! w_infty = rho_infty*k_infty/mu_t_ed_infty

    ! Farfield turbulence intensity
    ! k_infty = (3/2)*SSG_LRRW_tu_infty^2*U_infty^2
    !real(rk), parameter :: SSG_LRRW_tu_infty    = 0.00039_rk 
    real(rk), parameter :: sst_tu_infty    = 0.0001_rk 


    ! Farfield eddy viscosity fraction
    ! mu_t_ed_infty = SSG_LRRW_f_mu_infty*mu_infty
    !real(rk), parameter :: SSG_LRRW_f_mu_infty  = 0.009_rk
    real(rk), parameter :: sst_f_mu_infty  = 0.1_rk

    real(rk), parameter :: sst_k_infty    = 1.0e-3_rk 
    real(rk), parameter :: sst_omega_infty    = 5.8e7_rk*sst_f_mu_infty*sst_k_infty
end module mod_sst

