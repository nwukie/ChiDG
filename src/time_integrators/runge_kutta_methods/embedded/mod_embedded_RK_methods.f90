module mod_define_embedded_RK_methods
! Module containing subroutines and functions defining the Runge-Kutta-Fehlberg (RKF45) embedded method and
! the Dormand-Prince (DP54) embedded method
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO
    implicit none
    
    
    
   !!
   !! Embedded Runge-Kutta (RK) methods use two different Runge-Kutta methods of different orders, p and p' (can be p + 1/p - 1)
   !! to estimate the local truncation error in the method being used to advance the solution and the estimate is then used to 
   !! adjust the time step size which leads to faster convergence
   !!
   !! In the RKF45 method, a 5th order RK method is used to estimate the error in a 4th order RK method
   !!
   !! In the DP54 method, a 4th order RK method is used to estimate the local error in a 5th order RK method
   !!
   !! Both the methods are defined by their Butcher tableaus (see the reference below for more information on the Butcher tableau
   !! and both the methods)
   !!
   !! Reference: Hairer, E., Norsett, S.P., Wanner, G., Solving Ordinary Differential Equations I, Vol. 8, 
   !!            Springer Series in Computational Mathematics, Springer-Verlag, 1991
   !!
   !! The subroutines in this module return the number of stages and the coefficients (as defined in the Butecher tableau) for
   !! both methods
   !!
contains



    !> Define the Runge-Kutta-Fehlberg method (RKF45 with 6 stages)
    !!
    !! @author Mayank Sharma
    !! @date   1/20/2017
    !!
    !! @param[out]  nstage - Number of stages in the RKF45 method
    !! @param[out]  p      - Order of the solution advancing Runge-Kutta method
    !! @param[out]  a      - Array of coefficients for calculating stagewise updates
    !! @param[out]  b      - Array of coefficients for summing stagewise updates to get final update (4th order)
    !! @param[out]  b_app  - Array of coefficients for summing stagewise updates for 5th order solution 
    !!
    !------------------------------------------------------------------------------------------------------------------
    subroutine RKF45(nstage,p,a,b,b_app)
        integer(ik),             intent(inout)      :: nstage, p
        real(rk),   allocatable, intent(inout)      :: a(:,:)
        real(rk),   allocatable, intent(inout)      :: b(:)
        real(rk),   allocatable, intent(inout)      :: b_app(:)

        !
        ! Set number of stages
        !
        nstage = 6

        p = 4

        if (allocated(a) .and. allocated(b) .and. allocated(b_app)) deallocate(a,b,b_app)
        allocate(a(nstage,nstage),b(nstage),b_app(nstage))

        !
        ! a is a lower triangular matrix
        !
        a = ZERO

        !
        ! Set the values of the remaining coefficients in a
        !
        a(2,1) = 1._rk/4._rk
        a(3,1) = 3._rk/32._rk; a(3,2) = 9._rk/32._rk
        a(4,1) = 1932._rk/2197._rk; a(4,2) = -7200._rk/2197._rk; a(4,3) = 7296._rk/2197._rk
        a(5,1) = 439._rk/216._rk; a(5,2) = -8._rk; a(5,3) = 3680._rk/513._rk; a(5,4) = -845._rk/4104._rk
        a(6,1) = -8._rk/27._rk; a(6,2) = 2._rk; a(6,3) = -3544._rk/2565._rk; a(6,4) = 1859._rk/4104._rk; a(6,5) = -11._rk/40._rk

        !
        ! Set the values of the coefficients in b
        !
        b = [25.0_rk/216.0_rk, 0._rk, 1408._rk/2565._rk, 2197._rk/4104._rk, -1._rk/5._rk, 0._rk]

        !
        ! Set the values of the coefficients in b_app
        !
        b_app = [16._rk/135._rk, 0._rk, 6656._rk/12825._rk, 28561._rk/56430._rk, -9._rk/50._rk, 2._rk/55._rk]

    end subroutine RKF45
    !******************************************************************************************************************
    
    
    
    
    !> Define the Dormand-Prince method (DP54 with 7 stages)
    !!
    !! @author Mayank Sharma
    !! @date   1/20/2017
    !!
    !! @param[out]  nstage - Number of stages in the DP54 method
    !! @param[out]  p      - Order of the solution advancing Runge-Kutta method
    !! @param[out]  a      - Array of coefficients for calculating stagewise updates
    !! @param[out]  b      - Array of coefficients for summing stagewise updates to get final update (5th order)
    !! @param[out]  b_app  - Array of coefficients for summing stagewise updates for 4th order method
    !!
    !------------------------------------------------------------------------------------------------------------------
    subroutine DP54(nstage,p,a,b,b_app)
        integer(ik),             intent(inout)      :: nstage, p
        real(rk),   allocatable, intent(inout)      :: a(:,:)
        real(rk),   allocatable, intent(inout)      :: b(:)
        real(rk),   allocatable, intent(inout)      :: b_app(:)

        !
        ! Set number of stages
        !
        nstage = 7

        p = 5

        if (allocated(a) .and. allocated(b) .and. allocated(b_app)) deallocate(a,b,b_app)
        allocate(a(nstage,nstage),b(nstage),b_app(nstage))

        !
        ! a is a lower triangular matrix
        !
        a = ZERO

        !
        ! Set the values of the remaining coefficients in a
        !
        a(2,1) = 1._rk/5._rk
        a(3,1) = 3._rk/40._rk; a(3,2) = 9._rk/40._rk
        a(4,1) = 44._rk/45._rk; a(4,2) = -56._rk/15._rk; a(4,3) = 32._rk/9._rk
        a(5,1) = 19372._rk/6561._rk; a(5,2) = -25360._rk/2187._rk; a(5,3) = 64448._rk/6561._rk; a(5,4) = -212._rk/729._rk
        a(6,1) = 9017._rk/3168._rk; a(6,2) = -355._rk/33._rk; a(6,3) = 46732._rk/5247._rk; a(6,4) = 49._rk/176._rk
        a(6,5) = -5103._rk/18656._rk
        a(7,1) = 35._rk/384._rk; a(7,2) = 0._rk; a(7,3) = 500._rk/1113._rk; a(7,4) = 125._rk/192._rk; a(7,5) = -2187._rk/6784._rk
        a(7,6) = 11._rk/84._rk

        !
        ! Set the values of the coefficient in b
        ! Note that DP54 uses the FSAL (first same as last) strategy where the new solution is added as the 7th stage
        !
        b = [35._rk/384._rk, 0._rk, 500._rk/1113._rk, 125._rk/192._rk, -2187._rk/6784._rk, 11._rk/84._rk, 0._rk]

        !
        ! Set the values of the coefficients in b_app
        !
        b_app = [5179._rk/57600._rk, 0._rk, 7571._rk/16695._rk, 393._rk/640._rk, -92097._rk/339200._rk, 187._rk/2100._rk, 1._rk/40._rk]

    end subroutine DP54
    !******************************************************************************************************************



    !> Select embedded Runge_kutta method based on input
    !!
    !! @author Mayank Sharma
    !! @date   1/23/2017
    !!
    !! @param[in]   time_scheme - String denoting the embedded method to be used from user input
    !1 @param[out]  nstage      - Number of stages in the selected embedded method
    !! @param[out]  p           - Order of the solution advancing Runge-Kutta method
    !! @param[out]  a           - Array of coefficients for calculating stagewise updates
    !! @param[out]  b           - Array of coefficients for summing stagewise updates of solutiona dvancing method
    !! @param[out]  b_app       - Array of coefficients for summing stagewise updates of error estimating method
    !!
    !------------------------------------------------------------------------------------------------------------------
    subroutine method_selector(time_scheme,nstage,p,a,b,b_app)
        character(len = :), allocatable, intent(in)     :: time_scheme
        integer(ik),                     intent(inout)  :: nstage, p
        real(rk),           allocatable, intent(inout)  :: a(:,:)
        real(rk),           allocatable, intent(inout)  :: b(:) 
        real(rk),           allocatable, intent(inout)  :: b_app(:)


        select case(time_scheme)
            
            case('Runge-Kutta-Fehlberg Method', 'RKF', 'RKF45')
                call RKF45(nstage,p,a,b,b_app)

            case('Dormand-Prince Method', 'DP54')
                call DP54(nstage,p,a,b,b_app)

        end select


    end subroutine method_selector
    !******************************************************************************************************************



























end module mod_define_embedded_RK_methods
