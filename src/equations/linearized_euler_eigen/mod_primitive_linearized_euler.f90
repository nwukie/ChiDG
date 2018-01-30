module mod_primitive_linearized_euler
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: HALF, ONE, TWO, ZERO, PI
    implicit none

    !
    ! Azimuthal orther
    !
    integer(ik) :: mod_m = 2    ! Could get set by boundary conditions
    integer(ik) :: mod_n = 1    ! Could get set by boundary conditions


    !
    ! Nondimensionalization
    !
    real(rk), parameter :: u_r    = 1._rk
    real(rk), parameter :: rho_r  = 1._rk
    real(rk), parameter :: l_r    = 1._rk


    real(rk), parameter :: gam       = 1.4_rk
    !real(rk), parameter :: omega = 1432.945684_rk * TWO * PI
    !real(rk), parameter :: omega     = 956._rk * TWO * PI
    !real(rk), parameter :: omega     = ZERO
    real(rk), parameter :: omega     = 3441.9548_rk
    real(rk), parameter :: thickness = 0.5_rk
    real(rk), parameter :: eps       = 2300_rk

    !
    ! Dimensional mean flow
    !
    real(rk), parameter :: rho_d = 1.2212179_rk
    real(rk), parameter :: u_d   = 0._rk
    real(rk), parameter :: v_d   = 0._rk
    real(rk), parameter :: w_d   = 103.2586448_rk
    real(rk), parameter :: p_d   = 103341.6659_rk


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
    ! A-Matrix
    !
    real(rk), parameter :: rho_1_rho = ubar
    real(rk), parameter :: rho_1_u   = rhobar
    real(rk), parameter :: rho_1_v   = ZERO
    real(rk), parameter :: rho_1_w   = ZERO
    real(rk), parameter :: rho_1_p   = ZERO

    real(rk), parameter :: u_1_rho   = ZERO
    real(rk), parameter :: u_1_u     = ubar
    real(rk), parameter :: u_1_v     = ZERO
    real(rk), parameter :: u_1_w     = ZERO
    real(rk), parameter :: u_1_p     = ONE/rhobar

    real(rk), parameter :: v_1_rho   = ZERO
    real(rk), parameter :: v_1_u     = ZERO
    real(rk), parameter :: v_1_v     = ubar
    real(rk), parameter :: v_1_w     = ZERO
    real(rk), parameter :: v_1_p     = ZERO

    real(rk), parameter :: w_1_rho   = ZERO
    real(rk), parameter :: w_1_u     = ZERO
    real(rk), parameter :: w_1_v     = ZERO
    real(rk), parameter :: w_1_w     = ubar
    real(rk), parameter :: w_1_p     = ZERO

    real(rk), parameter :: p_1_rho   = ZERO
    real(rk), parameter :: p_1_u     = gam*pbar
    real(rk), parameter :: p_1_v     = ZERO
    real(rk), parameter :: p_1_w     = ZERO
    real(rk), parameter :: p_1_p     = ubar






    !
    ! B-Matrix
    !
    real(rk), parameter :: rho_2_rho = vbar
    real(rk), parameter :: rho_2_u   = ZERO
    real(rk), parameter :: rho_2_v   = rhobar
    real(rk), parameter :: rho_2_w   = ZERO
    real(rk), parameter :: rho_2_p   = ZERO

    real(rk), parameter :: u_2_rho   = ZERO
    real(rk), parameter :: u_2_u     = vbar
    real(rk), parameter :: u_2_v     = ZERO
    real(rk), parameter :: u_2_w     = ZERO
    real(rk), parameter :: u_2_p     = ZERO

    real(rk), parameter :: v_2_rho   = ZERO
    real(rk), parameter :: v_2_u     = ZERO
    real(rk), parameter :: v_2_v     = vbar
    real(rk), parameter :: v_2_w     = ZERO
    real(rk), parameter :: v_2_p     = ONE/rhobar

    real(rk), parameter :: w_2_rho   = ZERO
    real(rk), parameter :: w_2_u     = ZERO
    real(rk), parameter :: w_2_v     = ZERO
    real(rk), parameter :: w_2_w     = vbar
    real(rk), parameter :: w_2_p     = ZERO

    real(rk), parameter :: p_2_rho   = ZERO
    real(rk), parameter :: p_2_u     = ZERO
    real(rk), parameter :: p_2_v     = gam*pbar
    real(rk), parameter :: p_2_w     = ZERO
    real(rk), parameter :: p_2_p     = vbar




    !
    ! C-Matrix
    !
    real(rk), parameter :: rho_3_rho = wbar
    real(rk), parameter :: rho_3_u   = ZERO
    real(rk), parameter :: rho_3_v   = ZERO
    real(rk), parameter :: rho_3_w   = rhobar
    real(rk), parameter :: rho_3_p   = ZERO

    real(rk), parameter :: u_3_rho   = ZERO
    real(rk), parameter :: u_3_u     = wbar
    real(rk), parameter :: u_3_v     = ZERO
    real(rk), parameter :: u_3_w     = ZERO
    real(rk), parameter :: u_3_p     = ZERO

    real(rk), parameter :: v_3_rho   = ZERO
    real(rk), parameter :: v_3_u     = ZERO
    real(rk), parameter :: v_3_v     = wbar
    real(rk), parameter :: v_3_w     = ZERO
    real(rk), parameter :: v_3_p     = ZERO

    real(rk), parameter :: w_3_rho   = ZERO
    real(rk), parameter :: w_3_u     = ZERO
    real(rk), parameter :: w_3_v     = ZERO
    real(rk), parameter :: w_3_w     = wbar
    real(rk), parameter :: w_3_p     = ONE/rhobar

    real(rk), parameter :: p_3_rho   = ZERO
    real(rk), parameter :: p_3_u     = ZERO
    real(rk), parameter :: p_3_v     = ZERO
    real(rk), parameter :: p_3_w     = gam*pbar
    real(rk), parameter :: p_3_p     = wbar



!    !
!    ! D-Matrix
!    !
!    real(rk), parameter :: s1_rho = ZERO
!    real(rk), parameter :: s1_u   = ZERO
!    real(rk), parameter :: s1_v   = ZERO
!    real(rk), parameter :: s1_w   = ZERO
!    real(rk), parameter :: s1_p   = ZERO
!
!    real(rk), parameter :: s2_rho = vbar*vbar/(rhobar*r)
!    real(rk), parameter :: s2_u   = wbar
!    real(rk), parameter :: s2_v   = ZERO
!    real(rk), parameter :: s2_w   = ZERO
!    real(rk), parameter :: s2_p   = ZERO
!
!    real(rk), parameter :: s3_rho = ZERO
!    real(rk), parameter :: s3_u   = ZERO
!    real(rk), parameter :: s3_v   = wbar
!    real(rk), parameter :: s3_w   = ZERO
!    real(rk), parameter :: s3_p   = ZERO
!
!    real(rk), parameter :: s4_rho = ZERO
!    real(rk), parameter :: s4_u   = ZERO
!    real(rk), parameter :: s4_v   = ZERO
!    real(rk), parameter :: s4_w   = wbar
!    real(rk), parameter :: s4_p   = ONE/rhobar
!
!    real(rk), parameter :: s5_rho = ZERO
!    real(rk), parameter :: s5_u   = ZERO
!    real(rk), parameter :: s5_v   = ZERO
!    real(rk), parameter :: s5_w   = gam*pbar
!    real(rk), parameter :: s5_p   = wbar








contains









end module mod_primitive_linearized_euler
