module mod_HB_matrices

!*****************************************************************************************
!*                                                                                       *
!* This module contains subroutines for calculating matrices necessary for the HB scheme *
!* inv_E      - inverse Fourier transform matrix                                         *
!* diff_inv_E - the partial derivative of the inverse Fourier transform matrix wrt time  *
!* E          - Fourier transform matrix                                                 *
!* D          - psuedo spectral operator                                                 *
!* For K number of HB frequencies, w and N number of HB time levels (N = 2K + 1), t      *
!* the sizes of all matrices considered here is NxN                                      *
!*                                                                                       * 
!* The general form of inv_E:                                                            *
!* ([1 sin(w_1*t_1) ... sin(w_K*t_1) cos(w_1*t_1) ... cos(w_K*t_1)],                     *
!*  [1 sin(w_1*t_2) ... sin(w_K*t_2) cos(w_1*t_2) ... cos(w_K*t_2)],                     *
!*                                  .                                                    *
!*                                  .                                                    *
!*  [1 sin(w_1*t_N) ... sin(w_K*t_N) cos(w_1*t_N) ... cos(w_K*t_N)])                     *
!*                                                                                       *
!* The general form of diff_inv_E:                                                       *
!* ([0 w_1*cos(w_1*t_1) ... w_K*cos(w_K*t_1) -w_1*sin(w_1*t_1) ... -w_K*sin(w_K*t_1)],   *
!*  [0 w_1*cos(w_1*t_2) ... w_K*cos(w_K*t_2) -w_1*sin(w_1*t_2) ... -w_K*sin(w_K*t_2)],   *
!*                                          .                                            *
!*                                          .                                            *
!*  [0 w_1*cos(w_1*t_N) ... w_K*cos(w_K*t_N) -w_1*sin(w_1*t_N) ... -w_K*sin(w_K*t_N)])   *
!*                                                                                       *
!*****************************************************************************************
    
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO,ONE
    use mod_inv,        only: inv



    implicit none

    

contains


    
    !> Calculate the inverse Fourier transform matrix, inv(E)
    !!
    !! @author Mayank Sharma
    !! @date   1/4/2017
    !!
    !! @param[in]   nfreq - Number of HB frequencies
    !! @param[in]   ntime - Number of HB time levels (= 2nfreq + 1)
    !! @param[in]   omega - Array containing HB frequency values (size nfreq)
    !! @param[in]   t     - Array containing HB time level values (size ntime)
    !!
    !------------------------------------------------------------------------------------
    subroutine calc_inv_E(nfreq,ntime,omega,t,inv_E)
        integer(ik),                    intent(in)          :: nfreq,ntime
        real(rk),dimension(:),          intent(in)          :: omega,t
        real(rk),dimension(ntime,ntime),intent(inout)       :: inv_E 

        integer                                             :: irow,icol

        
        
        !
        ! inv_E - inverse Fourier transform matrix
        !
        do irow = 1,ntime
            do icol = 1,ntime

                if (icol == 1) then

                    inv_E(irow,icol) = ONE

                else if (icol >= 2 .and. icol <= nfreq + 1) then

                    inv_E(irow,icol) = sin(omega(icol - 1)*t(irow))

                else 

                    inv_E(irow,icol) = cos(omega(icol - (nfreq + 1))*t(irow))

                end if

            end do  ! icol
        end do  ! irow


    end subroutine calc_inv_E
    !************************************************************************************



    !> Calculate the partial derivative of inv(E) wrt time
    !!
    !! @author Mayank Sharma
    !! @date   1/4/2017
    !!
    !! @param[in]   nfreq - Number of HB frequencies 
    !! @param[in]   ntime - Number of HB time levels (= 2nfreq + 1)
    !! @param[in]   omega - Array containing HB frequency values (size nfreq)
    !! @param[in]   t     - Array containing HB time levels (size ntime)
    !!
    !------------------------------------------------------------------------------------
    subroutine calc_diff_inv_E(nfreq,ntime,omega,t,diff_inv_E)
        integer(ik),                    intent(in)          :: nfreq,ntime
        real(rk),dimension(:),          intent(in)          :: omega,t
        real(rk),dimension(ntime,ntime),intent(inout)       :: diff_inv_E

        integer                                             :: irow,icol



        !
        ! diff_inv_E - partial derivative of inverse Fourier transform matrix wrt t
        !
        do irow = 1,ntime
            do icol = 1,ntime

                if (icol == 1) then

                    diff_inv_E(irow,icol) = ZERO

                else if (icol >= 2 .and. icol <= nfreq + 1) then

                    diff_inv_E(irow,icol) = omega(icol - 1)*cos(omega(icol - 1)*t(irow))

                else

                    diff_inv_E(irow,icol) = omega(icol - (nfreq + 1))*&
                                            sin(omega(icol - (nfreq + 1))*t(irow))

                end if

            end do  ! icol
        end do  ! irow


    end subroutine calc_diff_inv_E
    !************************************************************************************



    !> Calculate the pseudo spectral operator 
    !!
    !! @author Mayank Sharma
    !! @date   1/4/2017
    !!
    !! @param[in]   nfreq - Number of HB frequencies 
    !! @param[in]   ntime - Number of HB time levels (= 2nfreq + 1)
    !! @param[in]   omega - Array containing HB frequency values (size nfreq)
    !! @param[in]   t     - Array containing HB time levels (size ntime)
    !!
    !------------------------------------------------------------------------------------
    subroutine calc_pseudo_spectral_operator(nfreq,ntime,omega,t,D)
        integer(ik),                    intent(in)          :: nfreq,ntime
        real(rk),dimension(:),          intent(in)          :: omega,t
        real(rk),dimension(ntime,ntime),intent(inout)       :: D 

        real(rk),dimension(ntime,ntime)                     :: inv_E,diff_inv_E,E



        !
        ! D - pseudo spectral operator
        !
        call calc_inv_E(nfreq,ntime,omega,t,inv_E)

        call calc_diff_inv_E(nfreq,ntime,omega,t,diff_inv_E)

        E = inv(inv_E)

        D = matmul(diff_inv_E,E)


    end subroutine calc_pseudo_spectral_operator
    !************************************************************************************




















end module mod_HB_matrices
