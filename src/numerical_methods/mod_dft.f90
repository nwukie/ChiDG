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
    subroutine dft(inreal,inimag,outreal,outimag)
        type(AD_D), intent(in)                 :: inreal(:)
        type(AD_D), intent(in)                 :: inimag(:)
        type(AD_D), intent(inout), allocatable :: outreal(:)
        type(AD_D), intent(inout), allocatable :: outimag(:)

        real(rk)    :: theta, freq
        integer     :: n, imode, itheta

        ! Check odd-count input
        if (mod(size(inreal),2) == 0) call chidg_signal(FATAL,"dft: expecting odd-count of input values.")

        ! Initialize derivative arrays
        outreal = inreal
        outimag = inreal

        outreal = ZERO
        outimag = ZERO

        ! Loop over output modes
        n = size(inreal)
        do imode = 1,n
            do itheta = 1,n
                freq  = TWO*PI*real(imode-1,rk)
                theta = real(itheta-1,rk) / real(n,rk)
                outreal(imode) = outreal(imode)  +  inreal(itheta)*cos(freq*theta)  +  inimag(itheta)*sin(freq*theta)
                outimag(imode) = outimag(imode)  -  inreal(itheta)*sin(freq*theta)  +  inimag(itheta)*cos(freq*theta)
            end do !itheta
        end do !imode

        ! Normalize
        outreal = outreal / real(n,rk)
        outimag = outimag / real(n,rk)

    end subroutine dft
    !******************************************************************



    !>  Computes the partial discrete Fourier transform of an input data set.
    !!
    !!  Here, we wish to compute the dft of a signal defined from 0 to 2pi.
    !!  However, the incoming signal is defined only on a section of that 
    !!  range 0 to 2pi/T, where T is some integer indicating the number of
    !!  sections the original range has been broken into.
    !!
    !!  For mode numbers that are multiples of T, their dft can be computed 
    !!  as
    !!
    !!      F[k] = (T/N)sum[ f(i) e^(-j k 2pi n/N) ]
    !!
    !!  ADAPTED FROM:
    !!  www.nayuki.io/res/how-to-implement-the-discrete-fourier-transform/dft.py
    !!
    !!  LICENCE: 
    !!  Public Domain
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
    !!  @param[in]      factor      number of 'passages' in periodic signal 
    !!                              that the dft of the input signal will 
    !!                              have to be multiplied by.
    !!  @param[in]      modes       'm' in e^i m theta. Starting from 0.
    !!
    !-----------------------------------------------------------------
    subroutine pdft(inreal,inimag,outreal,outimag,factor,modes)
        type(AD_D),     intent(in)                  :: inreal(:)
        type(AD_D),     intent(in)                  :: inimag(:)
        type(AD_D),     intent(inout), allocatable  :: outreal(:)
        type(AD_D),     intent(inout), allocatable  :: outimag(:)
        integer(ik),    intent(in)                  :: factor
        integer(ik),    intent(in)                  :: modes(:)

        real(rk)    :: theta, mode
        integer     :: n, imode, itheta, ierr

        ! Check odd-count input
        if (mod(size(inreal),2) == 0) call chidg_signal(FATAL,"pdft: expecting odd-count of input values.")

        ! Initialize derivative arrays
        allocate(outreal(size(modes)), outimag(size(modes)), stat=ierr)
        if (ierr/=0) call AllocationError
        outreal = ZERO*inreal(1)
        outimag = ZERO*inreal(1)

        ! Loop over output modes
        n = size(inreal)
        do imode = 1,size(modes)
            mode = real(modes(imode),rk)
            do itheta = 1,n
                theta = TWO*PI*real(itheta-1,rk)/real(n*factor,rk)
                outreal(imode) = outreal(imode)  +  inreal(itheta)*cos(mode*theta)  +  inimag(itheta)*sin(mode*theta)
                outimag(imode) = outimag(imode)  -  inreal(itheta)*sin(mode*theta)  +  inimag(itheta)*cos(mode*theta)
            end do !itheta
        end do !imode

        ! Normalize
        outreal = real(factor,rk)*outreal/real(n*factor,rk)
        outimag = real(factor,rk)*outimag/real(n*factor,rk)

    end subroutine pdft
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
    !!  @param[in]  theta   location T normalized by the periodicity. T=t/P
    !!
    !!  For example, if the original dft was computed over a period P=2.5, then
    !!  say we want to evaluate the transform at a location t=0.7. The incoming
    !!  value for theta would be T=0.7/2.5
    !!
    !---------------------------------------------------------------------------
    subroutine idft_eval(inreal,inimag,theta,outreal,outimag) 
        type(AD_D),     intent(in)      :: inreal(:)
        type(AD_D),     intent(in)      :: inimag(:)
        real(rk),       intent(in)      :: theta(:)
        type(AD_D),     intent(inout)   :: outreal(:)
        type(AD_D),     intent(inout)   :: outimag(:)

        real(rk)    :: freq
        integer     :: itheta, imode, ierr

        ! Allocate and initialize derivatives
        outreal(:) = inreal(1)
        outimag(:) = inreal(1)
        outreal = ZERO
        outimag = ZERO

        do itheta = 1,size(theta)
            do imode = 1,size(inreal)
                if (imode <= ((size(inreal)-1)/2 + 1)) then
                    freq = TWO*PI*real(imode-1,rk)  ! positive frequencies
                else
                    freq = -TWO*PI*real(size(inreal)-imode+1,rk) ! negative frequencies
                end if
                outreal(itheta) = outreal(itheta) + inreal(imode)*cos(freq*theta(itheta)) - inimag(imode)*sin(freq*theta(itheta))
                outimag(itheta) = outimag(itheta) + inreal(imode)*sin(freq*theta(itheta)) + inimag(imode)*cos(freq*theta(itheta))
            end do
        end do 

    end subroutine idft_eval
    !******************************************************************




    !>  Computes the inverse discrete Fourier transform of a partial DFT
    !!  evaluated with pdft at user-specified locations, instead of the 
    !!  traditional equi-spaced grid.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/15/2018
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
    !---------------------------------------------------------------------------
    subroutine ipdft_eval(inreal,inimag,theta,factor,modes,outreal,outimag) 
        type(AD_D),     intent(in)      :: inreal(:)
        type(AD_D),     intent(in)      :: inimag(:)
        real(rk),       intent(in)      :: theta(:)
        integer(ik),    intent(in)      :: factor
        integer(ik),    intent(in)      :: modes(:)
        type(AD_D),     intent(inout)   :: outreal(:)
        type(AD_D),     intent(inout)   :: outimag(:)

        real(rk)    :: mode
        integer     :: itheta, imode, ierr

        ! Allocate and initialize derivatives
        outreal(:) = inreal(1)
        outimag(:) = inreal(1)
        outreal = ZERO
        outimag = ZERO

        do itheta = 1,size(theta)
            do imode = 1,size(inreal)
                mode = real(modes(imode),rk)
                outreal(itheta) = outreal(itheta) + inreal(imode)*cos(mode*TWO*PI*theta(itheta)) - inimag(imode)*sin(mode*TWO*PI*theta(itheta))
                outimag(itheta) = outimag(itheta) + inreal(imode)*sin(mode*TWO*PI*theta(itheta)) + inimag(imode)*cos(mode*TWO*PI*theta(itheta))
            end do !imode
        end do !itheta

    end subroutine ipdft_eval
    !******************************************************************








end module mod_dft
