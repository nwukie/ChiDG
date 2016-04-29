module mod_linearized_euler
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: HALF, ONE, TWO, ZERO
    implicit none

    !
    ! Nondimensionalization
    !
    real(rk), parameter :: u_r    = 1._rk
    real(rk), parameter :: rho_r  = 1._rk
    real(rk), parameter :: l_r    = 1._rk


    real(rk), parameter :: gam = 1.4_rk

    !
    ! Dimensional mean flow
    !
    real(rk), parameter :: rho_d  = 1.2_rk
    real(rk), parameter :: rhou_d = 0._rk
    real(rk), parameter :: rhov_d = 0._rk
    real(rk), parameter :: rhow_d = 0._rk
    real(rk), parameter :: rhoE_d = 260000.0_rk


    !
    ! Nondimensional mean flow
    !
    real(rk), parameter :: rho_c  = rho_d  * (l_r/(rho_r*u_r))
    real(rk), parameter :: rhou_c = rhou_d * (l_r/(rho_r*u_r*u_r))
    real(rk), parameter :: rhov_c = rhov_d * (l_r/(rho_r*u_r*u_r))
    real(rk), parameter :: rhow_c = rhow_d * (l_r/(rho_r*u_r*u_r))
    real(rk), parameter :: rhoE_c = rhoE_d * (l_r/(rho_r*u_r*u_r*u_r))


    !
    ! Mean primitives
    !
    real(rk), parameter :: rhobar = rho_c
    real(rk), parameter :: ubar   = rhou_c / rho_c
    real(rk), parameter :: vbar   = rhov_c / rho_c
    real(rk), parameter :: wbar   = rhow_c / rho_c
    real(rk), parameter :: pbar   = (gam - ONE) * (rhoE_c - HALF*( (rhou_c*rhou_c) + (rhov_c*rhov_c) + (rhow_c*rhow_c))/rho_c )
    real(rk), parameter :: Hbar   = (rhoE_c + pbar) / rho_c
    real(rk), parameter :: cbar   = sqrt(gam * pbar / rhobar)



    !
    ! Jacobians
    !
    real(rk), parameter :: dp_drho  = ((gam - ONE)/TWO) * ( ubar**TWO + vbar**TWO )
    real(rk), parameter :: dp_drhou = -(gam - ONE)*ubar
    real(rk), parameter :: dp_drhov = -(gam - ONE)*vbar
    real(rk), parameter :: dp_drhow = -(gam - ONE)*wbar
    real(rk), parameter :: dp_drhoE =  (gam - ONE)

   

    
    real(rk), parameter :: rho_x_rho  = 0._rk
    real(rk), parameter :: rho_x_rhou = 1._rk
    real(rk), parameter :: rho_x_rhov = 0._rk
    real(rk), parameter :: rho_x_rhow = 0._rk
    real(rk), parameter :: rho_x_rhoE = 0._rk

    real(rk), parameter :: rhou_x_rho  = -ubar**TWO + dp_drho
    real(rk), parameter :: rhou_x_rhou = TWO*ubar + dp_drhou
    real(rk), parameter :: rhou_x_rhov = dp_drhov
    real(rk), parameter :: rhou_x_rhow = dp_drhow
    real(rk), parameter :: rhou_x_rhoE = dp_drhoE

    real(rk), parameter :: rhov_x_rho  = -ubar * vbar
    real(rk), parameter :: rhov_x_rhou = vbar
    real(rk), parameter :: rhov_x_rhov = ubar
    real(rk), parameter :: rhov_x_rhow = ZERO
    real(rk), parameter :: rhov_x_rhoE = ZERO

    real(rk), parameter :: rhow_x_rho  = -ubar * wbar
    real(rk), parameter :: rhow_x_rhou = wbar
    real(rk), parameter :: rhow_x_rhov = ZERO
    real(rk), parameter :: rhow_x_rhow = ubar
    real(rk), parameter :: rhow_x_rhoE = ZERO

    real(rk), parameter :: rhoE_x_rho  = ubar*(dp_drho - Hbar)
    real(rk), parameter :: rhoE_x_rhou = Hbar + (ubar*dp_drhou)
    real(rk), parameter :: rhoE_x_rhov = ubar*dp_drhov
    real(rk), parameter :: rhoE_x_rhow = ubar*dp_drhow
    real(rk), parameter :: rhoE_x_rhoE = ubar*(ONE + dp_drhoE)




    real(rk), parameter :: rho_y_rho  = 0._rk
    real(rk), parameter :: rho_y_rhou = 0._rk
    real(rk), parameter :: rho_y_rhov = 1._rk
    real(rk), parameter :: rho_y_rhow = 0._rk
    real(rk), parameter :: rho_y_rhoE = 0._rk

    real(rk), parameter :: rhou_y_rho  = -ubar * vbar
    real(rk), parameter :: rhou_y_rhou = vbar
    real(rk), parameter :: rhou_y_rhov = ubar
    real(rk), parameter :: rhou_y_rhow = ZERO
    real(rk), parameter :: rhou_y_rhoE = ZERO

    real(rk), parameter :: rhov_y_rho  = -vbar**TWO + dp_drho
    real(rk), parameter :: rhov_y_rhou = dp_drhou
    real(rk), parameter :: rhov_y_rhov = TWO*vbar + dp_drhov
    real(rk), parameter :: rhov_y_rhow = dp_drhow
    real(rk), parameter :: rhov_y_rhoE = dp_drhoE

    real(rk), parameter :: rhow_y_rho  = -vbar*wbar
    real(rk), parameter :: rhow_y_rhou = ZERO
    real(rk), parameter :: rhow_y_rhov = wbar
    real(rk), parameter :: rhow_y_rhow = vbar
    real(rk), parameter :: rhow_y_rhoE = ZERO

    real(rk), parameter :: rhoE_y_rho  = vbar*(dp_drho - Hbar)
    real(rk), parameter :: rhoE_y_rhou = vbar*dp_drhou
    real(rk), parameter :: rhoE_y_rhov = Hbar + (vbar*dp_drhov)
    real(rk), parameter :: rhoE_y_rhow = vbar*dp_drhow
    real(rk), parameter :: rhoE_y_rhoE = vbar*(ONE + dp_drhoE)






    real(rk), parameter :: rho_z_rho  = 0._rk
    real(rk), parameter :: rho_z_rhou = 0._rk
    real(rk), parameter :: rho_z_rhov = 0._rk
    real(rk), parameter :: rho_z_rhow = 1._rk
    real(rk), parameter :: rho_z_rhoE = 0._rk

    real(rk), parameter :: rhou_z_rho  = -ubar*wbar
    real(rk), parameter :: rhou_z_rhou = wbar
    real(rk), parameter :: rhou_z_rhov = ZERO
    real(rk), parameter :: rhou_z_rhow = ubar
    real(rk), parameter :: rhou_z_rhoE = ZERO

    real(rk), parameter :: rhov_z_rho  = -vbar*wbar
    real(rk), parameter :: rhov_z_rhou = ZERO
    real(rk), parameter :: rhov_z_rhov = wbar
    real(rk), parameter :: rhov_z_rhow = vbar
    real(rk), parameter :: rhov_z_rhoE = ZERO

    real(rk), parameter :: rhow_z_rho  = -wbar**TWO + dp_drho
    real(rk), parameter :: rhow_z_rhou = dp_drhou
    real(rk), parameter :: rhow_z_rhov = dp_drhov
    real(rk), parameter :: rhow_z_rhow = TWO*wbar + dp_drhow
    real(rk), parameter :: rhow_z_rhoE = dp_drhoE

    real(rk), parameter :: rhoE_z_rho  = wbar*(dp_drho - Hbar)
    real(rk), parameter :: rhoE_z_rhou = wbar*dp_drhou
    real(rk), parameter :: rhoE_z_rhov = wbar*dp_drhov
    real(rk), parameter :: rhoE_z_rhow = wbar*dp_drhow + Hbar
    real(rk), parameter :: rhoE_z_rhoE = wbar*(ONE+dp_drhoE)




















contains









end module
