module mod_gridspace
#include <messenger.h>
    use mod_kinds,  only: rk, ik
    implicit none









contains



    !>  Compute a linear distribution of a given number of values between
    !!  two upper and lower values.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/14/2016
    !!
    !------------------------------------------------------------------------
    function linspace(xmin,xmax,nvals) result(vals)
        real(rk),       intent(in)  :: xmin
        real(rk),       intent(in)  :: xmax
        integer(ik),    intent(in)  :: nvals

        real(rk), dimension(nvals)  :: vals
        real(rk)                    :: interval, dval
        integer(ik)                 :: ival

        if ( nvals < 2 ) call chidg_signal(FATAL,"linspace: number of grid-space values is fewer than two.")

        !
        ! Compute linear interval
        !
        interval = xmax-xmin
        dval     = interval/real((nvals-1),rk)

        
        do ival = 1,nvals
            vals(ival) = xmin + real((ival-1),rk)*dval
        end do

    end function linspace
    !************************************************************************


end module mod_gridspace
