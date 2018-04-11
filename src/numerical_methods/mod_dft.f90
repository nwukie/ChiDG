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
    !!  NOTE: The normalization convention used here, (1/N), defines
    !!  the DC part of the transform X[0] to be the average.
    !!
    !!  The transform result is in 'standard' order:
    !!  [ DC, omega_1, omega_2,..., omega_n, -omega_n, ... -omega_2, -omega_1]
    !!
    !!  [DC|     positive frequencies      |       negative frequencies      ]
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/26/2018
    !!
    !!  @param[in]      input       periodic signal
    !!  @param[inout]   outreal     real part of DFT
    !!  @param[inout]   outimag     imag part of DFT
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
            do imode = 1,size(inreal)
                freq = TWO*PI*real(imode-1,rk)
                output(itheta) = output(itheta) + inreal(imode)*cos(freq*theta(itheta)) - inimag(imode)*sin(freq*theta(itheta))
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
    !!  @param[in]  theta   Array of locations where the idft will be evaluated 
    !!                      normalized by the period.
    !!  @param[in]  period  Period, which will be used to normalize the 'theta' 
    !!                      inputs with respect to the Fourier transform.
    !!
    !!  For example, if the original dft was computed over a period 2.5, then
    !!  say we want to evaluate the transform at a location 0.7. The incoming
    !!  value for theta would be 0.7/2.5
    !!
    !!  Incoming transform should be in 'standard' order:
    !!  [ DC, omega_1, omega_2,..., omega_n, -omega_n, ... -omega_2, -omega_1]
    !!
    !!  [DC|     positive frequencies      |       negative frequencies      ]
    !!
    !---------------------------------------------------------------------------
    subroutine idft_eval(inreal,inimag,theta,outreal,outimag) 
        type(AD_D),     intent(in)      :: inreal(:)
        type(AD_D),     intent(in)      :: inimag(:)
        real(rk),       intent(in)      :: theta(:)
        type(AD_D),     intent(inout)   :: outreal(:)
        type(AD_D),     intent(inout)   :: outimag(:)

        !type(AD_D), allocatable :: output(:)
        real(rk)                :: freq
        integer                 :: itheta, imode, ierr

        ! Check odd-count input
        if (mod(size(inreal),2) == 0) call chidg_signal(FATAL,"idft_eval: expecting odd-count of input values.")

        ! Allocate and initialize derivatives
        !allocate(outreal(size(theta)),outimag(size(theta)), stat=ierr)
        !if (ierr /= 0) call AllocationError
        outreal(:) = inreal(1)
        outimag(:) = inreal(1)
        outreal = ZERO
        outimag = ZERO

        do itheta = 1,size(theta)
            do imode = 1,size(inreal)
                freq = TWO*PI*real(imode-1,rk)
                outreal(itheta) = outreal(itheta) + inreal(imode)*cos(freq*theta(itheta)) - inimag(imode)*sin(freq*theta(itheta))
                outimag(itheta) = outimag(itheta) + inreal(imode)*sin(freq*theta(itheta)) + inimag(imode)*cos(freq*theta(itheta))
            end do
        end do 

    end subroutine idft_eval
    !******************************************************************



end module mod_dft
