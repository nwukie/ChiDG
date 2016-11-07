module mod_spalart_allmaras
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: ONE, TWO, THREE
    implicit none



    real(rk), parameter :: SA_c_b1 = 0.1355_rk
    real(rk), parameter :: SA_c_b2 = 0.622_rk
    real(rk), parameter :: SA_kappa = 0.41_rk
    real(rk), parameter :: SA_sigma = TWO/THREE
    real(rk), parameter :: SA_c_w1  = SA_c_b1/(SA_kappa*SA_kappa)  + (ONE + SA_c_b2)/SA_sigma
    real(rk), parameter :: SA_c_w2  = 0.3_rk
    real(rk), parameter :: SA_c_w3  = TWO
    real(rk), parameter :: SA_c_v1  = 7.1_rk
    real(rk), parameter :: SA_c_v2  = 0.7_rk
    real(rk), parameter :: SA_c_v3  = 0.9_rk
    real(rk), parameter :: SA_c_t3  = 1.2_rk
    real(rk), parameter :: SA_c_t4  = 0.5_rk
    real(rk), parameter :: SA_c_n1  = 16._rk
    real(rk), parameter :: SA_Pr_t  = 0.9_rk
    real(rk), parameter :: SA_rlim  = 5._rk
    real(rk), parameter :: SA_b     = 100._rk









end module mod_spalart_allmaras
