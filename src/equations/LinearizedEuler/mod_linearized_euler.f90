module mod_linearized_euler
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: HALF, ONE, TWO, ZERO
    implicit none


    real(rk), parameter :: gam = 1.4_rk

    !real(rk), parameter :: rho_c  = 1.1928790357326_rk
    real(rk), parameter :: rho_c  = 1._rk
    !real(rk), parameter :: rhou_c = 151.826008671035_rk
    !real(rk), parameter :: rhou_c = 50._rk
    real(rk), parameter :: rhou_c = 0.2_rk
    !real(rk), parameter :: rhou_c = 0.001_rk
    real(rk), parameter :: rhov_c = 0._rk
    real(rk), parameter :: rhow_c = 0._rk
    !real(rk), parameter :: rhoE_c = 259661.975866153_rk
    real(rk), parameter :: rhoE_c = 1._rk

    real(rk), parameter :: ubar = rhou_c / rho_c
    real(rk), parameter :: vbar = rhov_c / rho_c
    real(rk), parameter :: wbar = rhow_c / rho_c

    real(rk), parameter :: pbar = (gam - ONE) * (rhoE_c - HALF*( (rhou_c*rhou_c) + (rhov_c*rhov_c) + (rhow_c*rhow_c))/rho_c )
    real(rk), parameter :: Hbar = (rhoE_c + pbar) / rho_c

    real(rk), parameter :: dp_drho  = ((gam - ONE)/TWO) * ( ubar**TWO + vbar**TWO )
    real(rk), parameter :: dp_drhou = -(gam - ONE)*ubar
    real(rk), parameter :: dp_drhov = -(gam - ONE)*vbar
    real(rk), parameter :: dp_drhow = -(gam - ONE)*wbar
    real(rk), parameter :: dp_drhoE =  (gam - ONE)

   

    
    real(rk), parameter :: rho_x_rho  = 0._rk
    real(rk), parameter :: rho_x_rhou = 1._rk
    real(rk), parameter :: rho_x_rhov = 0._rk
    real(rk), parameter :: rho_x_rhoE = 0._rk

    real(rk), parameter :: rhou_x_rho  = -ubar**TWO + dp_drho
    real(rk), parameter :: rhou_x_rhou = TWO*ubar + dp_drhou
    real(rk), parameter :: rhou_x_rhov = dp_drhov
    real(rk), parameter :: rhou_x_rhoE = dp_drhoE

    real(rk), parameter :: rhov_x_rho  = -ubar * vbar
    real(rk), parameter :: rhov_x_rhou = vbar
    real(rk), parameter :: rhov_x_rhov = ubar
    real(rk), parameter :: rhov_x_rhoE = ZERO

    real(rk), parameter :: rhoE_x_rho  = ubar*(dp_drho - Hbar)
    real(rk), parameter :: rhoE_x_rhou = Hbar + (ubar*dp_drhou)
    real(rk), parameter :: rhoE_x_rhov = ubar*dp_drhov
    real(rk), parameter :: rhoE_x_rhoE = ubar*(ONE + dp_drhoE)




    real(rk), parameter :: rho_y_rho  = 0._rk
    real(rk), parameter :: rho_y_rhou = 0._rk
    real(rk), parameter :: rho_y_rhov = 1._rk
    real(rk), parameter :: rho_y_rhoE = 0._rk

    real(rk), parameter :: rhou_y_rho  = -ubar * vbar
    real(rk), parameter :: rhou_y_rhou = vbar
    real(rk), parameter :: rhou_y_rhov = ubar
    real(rk), parameter :: rhou_y_rhoE = ZERO

    real(rk), parameter :: rhov_y_rho  = -vbar**TWO + dp_drho
    real(rk), parameter :: rhov_y_rhou = dp_drhou
    real(rk), parameter :: rhov_y_rhov = TWO*vbar + dp_drhov
    real(rk), parameter :: rhov_y_rhoE = dp_drhoE

    real(rk), parameter :: rhoE_y_rho  = vbar*(dp_drho - Hbar)
    real(rk), parameter :: rhoE_y_rhou = vbar*dp_drhou
    real(rk), parameter :: rhoE_y_rhov = Hbar + (vbar*dp_drhov)
    real(rk), parameter :: rhoE_y_rhoE = vbar*(ONE + dp_drhoE)





contains









end module
