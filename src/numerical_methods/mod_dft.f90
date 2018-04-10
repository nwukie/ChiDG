!>  Discrete Fourier (Inverse)Transform utilities.
!!
!!  Procedures
!!  --------------------------------------
!!  dft         Compute the discrete Fourier transform of a signal.
!!  idft        Compute the inverse discrete Fourier transform of a frequency spectrum.
!!  idft_eval   Compute the inverse discrete Fourier transform evaluated at user-specified points.
!!
!!  @author Nathan A. Wukie
!!  @date   2/26/2018
!!
!----------------------------------------------------------------------
module mod_dft
#include <messenger.h>
    use mod_kinds,      only: ik, rk
    use mod_constants,  only: PI, ZERO, ONE, TWO
    use mod_gridspace,  only: linspace
    use DNAD_D
    implicit none


contains



    !>  Computes the discrete Fourier transform of an input data set.
    !!
    !!  ADAPTED FROM:
    !!  www.nayuki.io/res/how-to-implement-the-discrete-fourier-transform/dft.py
    !!
    !!  LICENCE: 
    !!  Public Domain
    !!
    !!  Computes X[i] = (1/N)sum[ x(i)e^(-j*2pi*n*k*/N) ]
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/26/2018
    !!
    !!  @param[in]      input       periodic signal
    !!  @param[inout]   outreal     real part of DFT
    !!  @param[inout]   outimag     imag part of DFT
    !!
    !!
    !-----------------------------------------------------------------
    subroutine dft(input,outreal,outimag)
        type(AD_D), intent(in)                 :: input(:)
        type(AD_D), intent(inout), allocatable :: outreal(:)
        type(AD_D), intent(inout), allocatable :: outimag(:)

        real(rk)    :: theta, freq
        integer     :: n, imode, itheta

        ! Check odd-count input
        if (mod(size(input),2) == 0) call chidg_signal(FATAL,"dft: expecting odd-count of input values.")

        ! Initialize derivative arrays
        outreal = input
        outimag = input

        outreal = ZERO
        outimag = ZERO

        ! Loop over output modes
        n = size(input)
        do imode = 1,n
            do itheta = 1,n
                freq  = TWO*PI*real(imode-1,rk)
                theta = real(itheta-1,rk) / real(n,rk)
                outreal(imode) = outreal(imode)  +  input(itheta)*cos(freq*theta)
                outimag(imode) = outimag(imode)  -  input(itheta)*sin(freq*theta)
            end do !itheta
        end do !imode

        ! Normalize
        outreal = outreal / real(n,rk)
        outimag = outimag / real(n,rk)

    end subroutine dft
    !******************************************************************




    !>  Computes the inverse discrete Fourier transform of an input data set.
    !!
    !!  Computes x[i] = sum[ X(i)e^(j*2pi*n*k*/N) ]
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/26/2018
    !!
    !-----------------------------------------------------------------
    function idft(inreal,inimag) result(output)
        type(AD_D),     intent(in)  :: inreal(:)
        type(AD_D),     intent(in)  :: inimag(:)

        type(AD_D), allocatable :: output(:)
        real(rk),   allocatable :: theta(:)
        real(rk)                :: freq
        integer                 :: itheta, imode, ierr

        ! Check odd-count input
        if (mod(size(inreal),2) == 0) call chidg_signal(FATAL,"idft: expecting odd-count of input values.")

        ! Allocate and initialize derivatives
        allocate(output(size(inreal)), stat=ierr)
        if (ierr /= 0) call AllocationError
        output(:) = inreal(1)
        output = ZERO

        ! Theta at equispaced nodes
        theta = linspace(0._rk, real((size(inreal)-1),rk)/real(size(inreal)),size(inreal))

        ! Loop over output points
        do itheta = 1,size(theta)
            do imode = 1,(size(inreal)-1)/2
                freq = TWO*PI*real(imode-1,rk)
                if ( imode > 1 ) then
                    output(itheta) = output(itheta) + TWO * (inreal(imode)*cos(freq*theta(itheta)) - inimag(imode)*sin(freq*theta(itheta)))
                else
                    output(itheta) = output(itheta) + inreal(imode)*cos(freq*theta(itheta)) - inimag(imode)*sin(freq*theta(itheta))
                end if
            end do !imode
        end do !itheta

    end function idft
    !******************************************************************





    !>  Computes the inverse discrete Fourier transform evaluated
    !!  at user-specified locations, instead of the traditional
    !!  equi-spaced grid.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/26/2018
    !!
    !!  @param[in]  inreal  Real part of DFT. (could be directly from 'dft')
    !!  @param[in]  inimag  Imag part of DFT. (could be directly from 'dft')
    !!  @param[in]  theta   Array of locations where the idft will be evaluated.
    !!  @param[in]  period  Period, which will be used to normalize the 'theta' 
    !!                      inputs with respect to the Fourier transform.
    !!
    !---------------------------------------------------------------------------
    function idft_eval(inreal,inimag,theta,period) result(output)
        type(AD_D),     intent(in)  :: inreal(:)
        type(AD_D),     intent(in)  :: inimag(:)
        real(rk),       intent(in)  :: theta(:)
        real(rk),       intent(in)  :: period

        type(AD_D), allocatable :: output(:)
        real(rk),   allocatable :: theta_eval(:)
        real(rk)                :: freq
        integer                 :: itheta, imode, ierr

        ! Check odd-count input
        if (mod(size(inreal),2) == 0) call chidg_signal(FATAL,"idft_eval: expecting odd-count of input values.")

        ! Allocate and initialize derivatives
        allocate(output(size(theta)), stat=ierr)
        if (ierr /= 0) call AllocationError
        output(:) = inreal(1)
        output = ZERO

        ! Theta at equispaced nodes
        theta_eval = theta/period

        ! Loop over output points
        do itheta = 1,size(theta)
            do imode = 1,(size(inreal)-1)/2
                freq = TWO*PI*real(imode-1,rk)
                if ( imode > 1 ) then
                    output(itheta) = output(itheta) + TWO * (inreal(imode)*cos(freq*theta_eval(itheta)) - inimag(imode)*sin(freq*theta_eval(itheta)))
                else
                    output(itheta) = output(itheta) + inreal(imode)*cos(freq*theta_eval(itheta)) - inimag(imode)*sin(freq*theta_eval(itheta))
                end if
            end do !imode
        end do !itheta

    end function idft_eval
    !******************************************************************



end module mod_dft
