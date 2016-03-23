module mod_primitive_linearized_euler
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
    real(rk), parameter :: rho_d = 1.2_rk
    real(rk), parameter :: u_d   = 0._rk
    real(rk), parameter :: v_d   = 0._rk
    real(rk), parameter :: w_d   = 0._rk
    real(rk), parameter :: p_d   = 110000.0_rk


!    !
!    ! Nondimensional mean flow
!    !
!    real(rk), parameter :: rho_c  = rho_d  * (l_r/(rho_r*u_r))
!    real(rk), parameter :: rhou_c = rhou_d * (l_r/(rho_r*u_r*u_r))
!    real(rk), parameter :: rhov_c = rhov_d * (l_r/(rho_r*u_r*u_r))
!    real(rk), parameter :: rhow_c = rhow_d * (l_r/(rho_r*u_r*u_r))
!    real(rk), parameter :: rhoE_c = rhoE_d * (l_r/(rho_r*u_r*u_r*u_r))


    !
    ! Mean primitives
    !
    real(rk), parameter :: rhobar = rho_d
    real(rk), parameter :: ubar   = u_d
    real(rk), parameter :: vbar   = v_d
    real(rk), parameter :: wbar   = w_d
    real(rk), parameter :: pbar   = p_d
    real(rk), parameter :: cbar   = sqrt(gam * pbar / rhobar)
    real(rk), parameter :: Mbar   = ubar/cbar



    !
    ! X-Matrix
    !
    real(rk), parameter :: rho_x_rho = ubar
    real(rk), parameter :: rho_x_u   = rhobar
    real(rk), parameter :: rho_x_v   = ZERO
    real(rk), parameter :: rho_x_w   = ZERO
    real(rk), parameter :: rho_x_p   = ZERO

    real(rk), parameter :: u_x_rho   = ZERO
    real(rk), parameter :: u_x_u     = ubar
    real(rk), parameter :: u_x_v     = ZERO
    real(rk), parameter :: u_x_w     = ZERO
    real(rk), parameter :: u_x_p     = ONE/rhobar

    real(rk), parameter :: v_x_rho   = ZERO
    real(rk), parameter :: v_x_u     = ZERO
    real(rk), parameter :: v_x_v     = ubar
    real(rk), parameter :: v_x_w     = ZERO
    real(rk), parameter :: v_x_p     = ZERO

    real(rk), parameter :: w_x_rho   = ZERO
    real(rk), parameter :: w_x_u     = ZERO
    real(rk), parameter :: w_x_v     = ZERO
    real(rk), parameter :: w_x_w     = ubar
    real(rk), parameter :: w_x_p     = ZERO

    real(rk), parameter :: p_x_rho   = ZERO
    real(rk), parameter :: p_x_u     = gam*pbar
    real(rk), parameter :: p_x_v     = ZERO
    real(rk), parameter :: p_x_w     = ZERO
    real(rk), parameter :: p_x_p     = ubar






    !
    ! Y-Matrix
    !
    real(rk), parameter :: rho_y_rho = vbar
    real(rk), parameter :: rho_y_u   = ZERO
    real(rk), parameter :: rho_y_v   = rhobar
    real(rk), parameter :: rho_y_w   = ZERO
    real(rk), parameter :: rho_y_p   = ZERO

    real(rk), parameter :: u_y_rho   = ZERO
    real(rk), parameter :: u_y_u     = vbar
    real(rk), parameter :: u_y_v     = ZERO
    real(rk), parameter :: u_y_w     = ZERO
    real(rk), parameter :: u_y_p     = ZERO

    real(rk), parameter :: v_y_rho   = ZERO
    real(rk), parameter :: v_y_u     = ZERO
    real(rk), parameter :: v_y_v     = vbar
    real(rk), parameter :: v_y_w     = ZERO
    real(rk), parameter :: v_y_p     = ONE/rhobar

    real(rk), parameter :: w_y_rho   = ZERO
    real(rk), parameter :: w_y_u     = ZERO
    real(rk), parameter :: w_y_v     = ZERO
    real(rk), parameter :: w_y_w     = vbar
    real(rk), parameter :: w_y_p     = ZERO

    real(rk), parameter :: p_y_rho   = ZERO
    real(rk), parameter :: p_y_u     = ZERO
    real(rk), parameter :: p_y_v     = gam*pbar
    real(rk), parameter :: p_y_w     = ZERO
    real(rk), parameter :: p_y_p     = vbar




    !
    ! Z-Matrix
    !
    real(rk), parameter :: rho_z_rho = wbar
    real(rk), parameter :: rho_z_u   = ZERO
    real(rk), parameter :: rho_z_v   = ZERO
    real(rk), parameter :: rho_z_w   = rhobar
    real(rk), parameter :: rho_z_p   = ZERO

    real(rk), parameter :: u_z_rho   = ZERO
    real(rk), parameter :: u_z_u     = wbar
    real(rk), parameter :: u_z_v     = ZERO
    real(rk), parameter :: u_z_w     = ZERO
    real(rk), parameter :: u_z_p     = ZERO

    real(rk), parameter :: v_z_rho   = ZERO
    real(rk), parameter :: v_z_u     = ZERO
    real(rk), parameter :: v_z_v     = wbar
    real(rk), parameter :: v_z_w     = ZERO
    real(rk), parameter :: v_z_p     = ZERO

    real(rk), parameter :: w_z_rho   = ZERO
    real(rk), parameter :: w_z_u     = ZERO
    real(rk), parameter :: w_z_v     = ZERO
    real(rk), parameter :: w_z_w     = wbar
    real(rk), parameter :: w_z_p     = ONE/rhobar

    real(rk), parameter :: p_z_rho   = ZERO
    real(rk), parameter :: p_z_u     = ZERO
    real(rk), parameter :: p_z_v     = ZERO
    real(rk), parameter :: p_z_w     = gam*pbar
    real(rk), parameter :: p_z_p     = wbar











contains









end module mod_primitive_linearized_euler
