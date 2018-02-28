module mod_dft
#include <messenger.h>
    use mod_kinds,      only: ik, rk
    use mod_constants,  only: PI, ZERO, ONE, TWO
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
    !!  Computes X[i] = SUM[ x(i)e^(-i*2pi*n*k*/N) ]
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/12/2016
    !!
    !!
    !-----------------------------------------------------------------
    subroutine dft(input,outreal,outimag)
        type(AD_D), intent(in)                 :: input(:)
        type(AD_D), intent(inout), allocatable :: outreal(:)
        type(AD_D), intent(inout), allocatable :: outimag(:)

        type(AD_D) :: sreal, simag, theta, freq
        integer    :: n, imode, itheta


        !
        ! Initialize derivative arrays
        !
        outreal = input
        outimag = input
        sreal   = input(1)
        simag   = input(1)
        theta   = input(1)

        outreal = ZERO
        outimag = ZERO
        sreal   = ZERO
        simag   = ZERO
        theta   = ZERO



        !
        ! Loop over output modes
        !
        n = size(input)
        do imode = 1,n

            !
            ! Get contribution from input points
            !
            sreal = ZERO
            simag = ZERO
            do itheta = 1,n

                freq  = TWO*PI*real(imode-1,rk)
                theta = real(itheta-1,rk) / real(n,rk)

                sreal = sreal + input(itheta) * cos(-freq*theta)
                simag = simag + input(itheta) * sin(-freq*theta)

            end do !itheta

            outreal(imode) = sreal
            outimag(imode) = simag

        end do !imode

        outreal = outreal / real(n,rk)
        outimag = outimag / real(n,rk)

    end subroutine dft
    !******************************************************************









    !>  Computes the inverse discrete Fourier transform of an input data set.
    !!
    !!  Adapted from Public Domain code:
    !!  www.nayuki.io/res/how-to-implement-the-discrete-fourier-transform/dft.py
    !!
    !!  Computes x[i] = (1/N) SUM[ X(i)e^( i*2pi*n*k*/N) ]
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/12/2016
    !!
    !!
    !-----------------------------------------------------------------
    function idft(inreal,inimag,theta) result(output)
        type(AD_D),     intent(in)  :: inreal(:)
        type(AD_D),     intent(in)  :: inimag(:)
        real(rk),       intent(in)  :: theta(:)

        type(AD_D),   allocatable   :: output(:)
        type(AD_D)                  :: s

        real(rk)    :: freq
        integer     :: itheta, imode, ierr


        !
        ! Allocate and initialize derivatives
        !
        allocate(output(size(theta)), stat=ierr)
        if (ierr /= 0) call AllocationError
        output(:) = inreal(1)
        s         = inreal(1)

        output = ZERO
        s      = ZERO


        !
        ! Loop over output points
        !
        do itheta = 1,size(theta)
            output(itheta) = ZERO
            s              = ZERO
            do imode = 1,size(inreal)

                ! Compute frequency and location
                freq = TWO*PI*real(imode-1,rk)

                if ( imode > 1 ) then
                    s = s + TWO * (inreal(imode)*cos(freq*theta(itheta)) - inimag(imode)*sin(freq*theta(itheta)))
                else
                    s = s + inreal(imode)*cos(freq*theta(itheta)) - inimag(imode)*sin(freq*theta(itheta))
                end if
                
            end do !imode
            output(itheta) = s
        end do !itheta


    end function idft
    !******************************************************************






end module mod_dft
