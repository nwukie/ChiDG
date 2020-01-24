module mod_HB_matrices
#include<messenger.h>

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
        integer(ik),            intent(in)          :: nfreq,ntime
        real(rk),dimension(:),  intent(in)          :: omega,t
        real(rk),dimension(:,:),intent(inout)       :: inv_E 

        integer                                     :: irow,icol

        
        
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
        integer(ik),            intent(in)          :: nfreq,ntime
        real(rk),dimension(:),  intent(in)          :: omega,t
        real(rk),dimension(:,:),intent(inout)       :: diff_inv_E

        integer                                     :: irow,icol



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

                    diff_inv_E(irow,icol) = -omega(icol - (nfreq + 1))*&
                                            sin(omega(icol - (nfreq + 1))*t(irow))

                end if

            end do  ! icol
        end do  ! irow


    end subroutine calc_diff_inv_E
    !************************************************************************************



    !>  Calculate the Fourier transform matrix
    !!
    !!  @author Mayank Sharma
    !!  @date   3/18/2017
    !!
    !------------------------------------------------------------------------------------
    subroutine calc_E(nfreq,ntime,omega,t,E)
        integer(ik),    intent(in)      :: nfreq, ntime
        real(rk),       intent(in)      :: omega(:), t(:)
        real(rk),       intent(inout)   :: E(:,:)

        real(rk),   allocatable        :: inv_E(:,:)
        integer(ik)                    :: ierr


        !
        ! Compute the inverse Fourier transform matrix
        !
        if (allocated(inv_E)) deallocate(inv_E)
        allocate(inv_E(ntime,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError

        call calc_inv_E(nfreq,ntime,omega,t,inv_E)


        !
        ! E - Fourier transform matrix
        !
        E = inv(inv_E)


    end subroutine calc_E
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
    !subroutine calc_pseudo_spectral_operator(nfreq,ntime,omega,t,D)
    function calc_pseudo_spectral_operator(omega,t) result(D)
        real(rk),                   intent(in)      :: omega(:)
        real(rk),                   intent(in)      :: t(:)

        real(rk),       allocatable :: inv_E(:,:),diff_inv_E(:,:), E(:,:), D(:,:)
        character(:),   allocatable :: user_msg, dev_msg
        integer(ik)                 :: i,ierr, ntime, nfreq

        ntime = size(t)
        nfreq = size(omega)

        if (allocated(inv_E) .and. allocated(diff_inv_E) .and. allocated(E) &
            .and. allocated(D)) deallocate(inv_E,diff_inv_E,E,D)
        allocate(inv_E(ntime,ntime),diff_inv_E(ntime,ntime), &
                 E(ntime,ntime),D(ntime,ntime), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! D - pseudo spectral operator
        !
        call calc_inv_E(nfreq,ntime,omega,t,inv_E)
        call calc_diff_inv_E(nfreq,ntime,omega,t,diff_inv_E)
    
        user_msg = 'The size of an array being inverted here is ZERO. Check &
                    the frequencies you have entered in chidg.nml'
        dev_msg  = 'The inv subroutine expects a valid non-zero integer matrix size. Check &
                    in time_manager.f90'

        if (size(inv_E,1) == 0) call chidg_signal_two(FATAL, user_msg, size(inv_E,1), dev_msg)
        
        E = inv(inv_E)
        D = matmul(diff_inv_E,E)

    end function calc_pseudo_spectral_operator
    !************************************************************************************




















end module mod_HB_matrices
