module mod_interpolation
#include <messenger.h>
    use DNAD_D
    implicit none





contains


    !>  General interpolation routine. Based on given interpolation_type, select appropriate
    !!  interpolation routine(linear, quadratic, etc.).
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/2/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    !subroutine interpolate(interpolation_type, xarray, yarray, xnew, ynew)
    !    character(*),   intent(in)      :: interpolation_type
    !    real(rk),       intent(in)      :: xarray(:)
    !    real(rk),       intent(in)      :: yarray(:)
    !    real(rk),       intent(in)      :: xnew
    !    real(rk),       intent(inout)   :: ynew
    function interpolate(interpolation_type, xarray, yarray, xnew) result(ynew)
        character(*),   intent(in)      :: interpolation_type
        real(rk),       intent(in)      :: xarray(:)
        real(rk),       intent(in)      :: yarray(:)
        real(rk),       intent(in)      :: xnew

        real(rk) :: ynew


        if ( interpolation_type == 'linear' ) then
            !call interpolate_linear(xarray,yarray,xnew,ynew)
            ynew = interpolate_linear(xarray,yarray,xnew)

        else
            call chidg_signal_one(FATAL, "interpolate: invalid interpolation_type.", interpolation_type)

        end if



    end function interpolate
    !***************************************************************************************







    !>  Linearly interpolate a new value from a list of given x,y pairs.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/2/2016
    !!
    !---------------------------------------------------------------------------------------
    function interpolate_linear(xarray,yarray,xnew) result(ynew)
        real(rk),   intent(in)      :: xarray(:)
        real(rk),   intent(in)      :: yarray(:)
        real(rk),   intent(in)      :: xnew

        integer     :: ival, lindex, hindex
        real(rk)    :: xlow_current, xhigh_current, ynew
        logical     :: out_of_bounds


        !
        ! Test xnew is bounded by xarray min/max
        !
        out_of_bounds = ( xnew < minval(xarray) .or. xnew > maxval(xarray) )
        if ( out_of_bounds ) call chidg_signal_one(FATAL,"interpolate_linear: new x-value is not bounded by the min/max of given x-values.", xnew)


        !
        ! Find next lowest x-value
        !
        lindex          = minloc(xarray,dim=1)
        xlow_current    = minval(xarray)
        do ival = 1,size(xarray)

            if ( xarray(ival) < xnew  .and. xarray(ival) > xlow_current ) then
                xlow_current = xarray(ival)
                lindex       = ival
            end if

        end do



        !
        ! Find next highest x-value
        !
        hindex          = maxloc(xarray,dim=1)
        xhigh_current   = maxval(xarray)
        do ival = 1,size(xarray)

            if ( xarray(ival) > xnew  .and. xarray(ival) < xhigh_current ) then
                xhigh_current = xarray(ival)
                hindex        = ival
            end if

        end do

        
        !
        ! Compute linear interpolation between next highest/lowest values.
        !
        ynew = yarray(lindex) + ( yarray(hindex) - yarray(lindex) ) * ( xnew - xarray(lindex) )/( xarray(hindex) - xarray(lindex) )

    end function interpolate_linear
    !***************************************************************************************



    !>  Linearly interpolate a new value from a list of given x,y pairs.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/2/2016
    !!
    !---------------------------------------------------------------------------------------
    function interpolate_linear_ad(xarray,yarray,xnew) result(ynew)
        real(rk),   intent(in)      :: xarray(:)
        type(AD_D), intent(in)      :: yarray(:)
        real(rk),   intent(in)      :: xnew

        integer     :: ival, lindex, hindex
        real(rk)    :: xlow_current, xhigh_current
        logical     :: out_of_bounds
        type(AD_D)  :: ynew


        !
        ! Test xnew is bounded by xarray min/max
        !
        out_of_bounds = ( xnew < minval(xarray) .or. xnew > maxval(xarray) )
        if ( out_of_bounds ) call chidg_signal_one(FATAL,"interpolate_linear: new x-value is not bounded by the min/max of given x-values.", xnew)


        !
        ! Find next lowest x-value
        !
        lindex          = minloc(xarray,dim=1)
        xlow_current    = minval(xarray)
        do ival = 1,size(xarray)

            if ( xarray(ival) < xnew  .and. xarray(ival) > xlow_current ) then
                xlow_current = xarray(ival)
                lindex       = ival
            end if

        end do



        !
        ! Find next highest x-value
        !
        hindex          = maxloc(xarray,dim=1)
        xhigh_current   = maxval(xarray)
        do ival = 1,size(xarray)

            if ( xarray(ival) > xnew  .and. xarray(ival) < xhigh_current ) then
                xhigh_current = xarray(ival)
                hindex        = ival
            end if

        end do

        
        !
        ! Compute linear interpolation between next highest/lowest values.
        !
        ynew = yarray(lindex) + ( yarray(hindex) - yarray(lindex) ) * ( xnew - xarray(lindex) )/( xarray(hindex) - xarray(lindex) )

    end function interpolate_linear_ad
    !***************************************************************************************







end module mod_interpolation
