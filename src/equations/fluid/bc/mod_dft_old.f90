module mod_dft
#include <messenger.h>
    use mod_kinds,      only: ik, rk
    use mod_constants,  only: PI, ZERO, ONE, TWO, DIR_2, &
                              XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use type_domain,    only: domain_t
    use type_point,     only: point_t
    use DNAD_D
    implicit none







contains




    !>  Compute the points for a discrete Fourier transform across a boundary.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/17/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------
    function compute_dft_points(domain,elems,faces,periodicity) result(points)
        type(domain_t), intent(in)   :: domain
        integer(ik),    intent(in)   :: elems(:)
        integer(ik),    intent(in)   :: faces(:)
        real(rk),       intent(in)   :: periodicity

        type(point_t),  allocatable  :: points(:)
        integer(ik)                  :: nmodes, npoints, ierr, ielem_bc, ielem, var, mode, ipnt, min_y_element, elem_min, iface
        integer(ik)                  :: itime
        integer(ik)                  :: min_y_loc
        real(rk)                     :: xi, eta, zeta, min_y, dy, xloc, zloc

        real(rk),       allocatable  :: mean_y_coordinates(:)
        integer(ik),    allocatable  :: mean_y_elements(:)
        real(rk),       dimension(2) :: ylocs, xlocs, zlocs

        !
        ! TODO: fix hard coded face orientation
        !
        iface = faces(1)


        !
        ! Number of Fourier harmonics to compute
        !
        nmodes = 10


        !
        ! Compute number of points in dft
        !
        npoints = 1 + 2*nmodes


        !
        ! Compute distance between points
        !
        dy = periodicity/real((npoints), rk)


        !
        ! Allocate points
        !
        allocate(points(npoints), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Allocate storage for mean coordinates
        !
        allocate(mean_y_coordinates(size(elems)), &
                 mean_y_elements(   size(elems)), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Coordinate min/max on the boundary in direction of periodicity.
        !
        do ielem_bc = 1,size(elems)

            !
            ! Get element index in block domain
            !
            ielem = elems(ielem_bc)


            !
            ! Get mean y-coordinate from element. y-coordinate variable = 2. Mode = 1
            ! 
            var  = DIR_2    ! y-coordinate index
            mode = 1        ! Coordinate average
            mean_y_coordinates(ielem_bc) = domain%elems(ielem)%coords%getterm(var,mode,itime)
            mean_y_elements(ielem_bc)    = ielem

        end do ! ielem


        !
        ! Find element with minimum mean y-coordinate
        !
        min_y_element = minloc(mean_y_coordinates, DIM=1)


        elem_min = mean_y_elements(min_y_element)

        !
        ! For the element with the minimum mean y-coordinate, find minimum y-coordinate
        ! on the boundary face. So, for the time being, XI_MIN, ETA_MIN, ZETA=0
        !
!        !xi   = -ONE ! INLET
!        xi   = -ONE  ! OUTLET
!        eta  = -ONE
!        zeta = domain%elems(elem_min)%coords%getterm(3,1) ! get mean z-coordinate
!
!        min_y = domain%elems(elem_min)%y(xi,eta,zeta) + 0.0000001_rk


        if ( (iface == XI_MIN) ) then
            ylocs(1) = domain%elems(elem_min)%y(-ONE, ONE,ZERO)
            ylocs(2) = domain%elems(elem_min)%y(-ONE,-ONE,ZERO)

            xlocs(1) = domain%elems(elem_min)%x(-ONE, ONE,ZERO)
            xlocs(2) = domain%elems(elem_min)%x(-ONE,-ONE,ZERO)

            zlocs(1) = domain%elems(elem_min)%z(-ONE, ONE,ZERO)
            zlocs(2) = domain%elems(elem_min)%z(-ONE,-ONE,ZERO)

            min_y_loc = minloc(ylocs,1)
            min_y     = ylocs(min_y_loc) + 0.000001_rk
            xloc      = xlocs(min_y_loc)
            zloc      = zlocs(min_y_loc)

        else if ( (iface == XI_MAX) ) then
            ylocs(1) = domain%elems(elem_min)%y(ONE, ONE,ZERO)
            ylocs(2) = domain%elems(elem_min)%y(ONE,-ONE,ZERO)

            xlocs(1) = domain%elems(elem_min)%x(ONE, ONE,ZERO)
            xlocs(2) = domain%elems(elem_min)%x(ONE,-ONE,ZERO)

            zlocs(1) = domain%elems(elem_min)%z(ONE, ONE,ZERO)
            zlocs(2) = domain%elems(elem_min)%z(ONE,-ONE,ZERO)


            min_y_loc = minloc(ylocs,1)
            min_y     = ylocs(min_y_loc) + 0.000001_rk
            xloc      = xlocs(min_y_loc)
            zloc      = zlocs(min_y_loc)

        else if ( (iface == ETA_MIN) ) then
            ylocs(1) = domain%elems(elem_min)%y( ONE,-ONE,ZERO)
            ylocs(2) = domain%elems(elem_min)%y(-ONE,-ONE,ZERO)

            xlocs(1) = domain%elems(elem_min)%x( ONE,-ONE,ZERO)
            xlocs(2) = domain%elems(elem_min)%x(-ONE,-ONE,ZERO)

            zlocs(1) = domain%elems(elem_min)%z( ONE,-ONE,ZERO)
            zlocs(2) = domain%elems(elem_min)%z(-ONE,-ONE,ZERO)

            min_y_loc = minloc(ylocs,1)
            min_y     = ylocs(min_y_loc) + 0.000001_rk
            xloc      = xlocs(min_y_loc)
            zloc      = zlocs(min_y_loc)

        else if ( (iface == ETA_MAX) ) then
            ylocs(1) = domain%elems(elem_min)%y( ONE,ONE,ZERO)
            ylocs(2) = domain%elems(elem_min)%y(-ONE,ONE,ZERO)

            xlocs(1) = domain%elems(elem_min)%x( ONE,ONE,ZERO)
            xlocs(2) = domain%elems(elem_min)%x(-ONE,ONE,ZERO)

            zlocs(1) = domain%elems(elem_min)%z( ONE,ONE,ZERO)
            zlocs(2) = domain%elems(elem_min)%z(-ONE,ONE,ZERO)

            min_y_loc = minloc(ylocs,1)
            min_y     = ylocs(min_y_loc) + 0.000001_rk
            xloc      = xlocs(min_y_loc)
            zloc      = zlocs(min_y_loc)

        else if ( (iface == ZETA_MIN) .or. (iface == ZETA_MAX) ) then

            call chidg_signal(FATAL,"compute_dft_ponts: assumes boundary is on a XI or ETA face.")

        end if



        !
        ! Create array of points from (min_y) to (min_y + periodicity)
        !
        !call points(1)%set(xloc, min_y, ZERO)
        call points(1)%set(xloc, min_y, zloc)
        do ipnt = 2,npoints


            !
            ! Copy previous point
            !
            points(ipnt) = points(ipnt - 1)

            !
            ! Add dy
            !
            call points(ipnt)%add_y(dy)

        end do ! ipnt


    end function compute_dft_points
    !**********************************************************************************










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

        type(AD_D) :: sreal, simag, angle
        integer    :: n, k, t


        !
        ! Initialize derivative arrays
        !
        outreal = input
        outimag = input
        sreal   = input(1)
        simag   = input(1)
        angle   = input(1)

        outreal = ZERO
        outimag = ZERO
        sreal   = ZERO
        simag   = ZERO
        angle   = ZERO



        !
        ! Loop over output modes
        !
        n = size(input)
        !do k = 1,n+1
        do k = 1,n

            !
            ! Get contribution from input points
            !
            sreal = ZERO
            simag = ZERO
            do t = 1,n

                angle = TWO * PI * real(t-1,rk) * real(k-1,rk) / real(n,rk)

                sreal = sreal + input(t) * cos(-angle)
                simag = simag + input(t) * sin(-angle)

            end do

            outreal(k) = sreal
            outimag(k) = simag

        end do ! k


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
    subroutine idft(inreal,inimag,output)
        type(AD_D),   intent(in)                   :: inreal(:)
        type(AD_D),   intent(in)                   :: inimag(:)
        type(AD_D),   intent(inout), allocatable   :: output(:)

        type(AD_D)    :: s
        real(rk)       :: angle
        integer       :: k, t, nmodes, npoints


        nmodes   = size(inreal)
        npoints  = nmodes


        !
        ! Initialize derivative arrays
        !
        output = inreal
        s      = inreal(1)

        output = ZERO
        s      = ZERO


        !
        ! Loop over output points
        !
        do k = 1,npoints

            !
            ! Loop over input modes
            !
            output(k) = ZERO
            s         = ZERO
            do t = 1,nmodes

                angle = TWO * PI * real(t-1,rk) * real(k-1,rk) / real(nmodes,rk)

                s = s + inreal(t)*cos(angle) - inimag(t)*sin(angle)
                
            end do

            output(k) = s

        end do


        !
        ! Normalization
        !
        !output = output / real(nmodes,8)

    end subroutine idft
    !******************************************************************











    !>  Computes the inverse discrete Fourier transform of an input data at 
    !!  user specified angles.
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
    subroutine idft_angle(inreal, inimag, angles, output)
        type(AD_D),       intent(in)      :: inreal(:)
        type(AD_D),       intent(in)      :: inimag(:)
        type(AD_D),       intent(in)      :: angles(:)
        type(AD_D),       intent(inout), allocatable   :: output(:)

        type(AD_D)            :: sreal, simag, s, freq, mag, phase
        integer               :: iang, imode, nmodes, nangles, quadrant


        nangles  = size(angles)
        nmodes   = size(inreal)


        !
        ! Initialize derivative arrays
        !
        output = angles
        sreal  = inreal(1)
        simag  = inreal(1)
        s      = inreal(1)
        freq   = inreal(1)

        output = ZERO
        sreal  = ZERO
        simag  = ZERO
        s      = ZERO
        freq   = ZERO




        !
        ! Loop over angles
        !
        do iang = 1,nangles

            !
            ! Accumulate contribution from each mode.
            !
            output(iang) = ZERO
            s            = ZERO
            do imode = 1,(nmodes-1)/2

                
                !
                ! Compute current frequency
                !
                freq  = TWO * PI * real(imode-1,rk)

!                !
!                ! Compute magnitude and phase from the complex mode
!                !
!                mag   = sqrt( inreal(imode)**TWO + inimag(imode)**TWO )
!                phase = asin( inimag(imode) / mag )
!
!                !
!                ! Check quadrant
!                !
!                if      ( (inreal(imode) >= ZERO) .and. (inimag(imode) >= ZERO) ) then
!                    quadrant = 1
!                else if ( (inreal(imode) < ZERO)  .and. (inimag(imode) >= ZERO) ) then
!                    quadrant = 2
!                else if ( (inreal(imode) < ZERO)  .and. (inimag(imode) < ZERO) ) then
!                    quadrant = 3
!                else if ( (inreal(imode) >= ZERO)  .and. (inimag(imode) < ZERO) ) then
!                    quadrant = 4
!                end if
!
!
!
!
!                !
!                ! Account for quadrant
!                !
!                if      ( quadrant == 1 ) then
!                    phase = phase
!                else if ( quadrant == 2) then
!                    phase = PI - phase
!                else if ( quadrant == 3) then
!                    phase = PI - phase
!                else if ( quadrant == 4) then
!                    phase = phase
!                end if
!
!
!
!
!                if ( imode > 1 ) then
!                    s = s +  TWO * mag*cos( freq * angles(iang) + phase )
!                else
!                    s = s +  mag*cos( freq * angles(iang) + phase )
!                end if


                if ( imode > 1 ) then
                    s = s + TWO * (inreal(imode)*cos(freq*angles(iang)) - inimag(imode)*sin(freq*angles(iang)))
                else
                    s = s + inreal(imode)*cos(freq*angles(iang)) - inimag(imode)*sin(freq*angles(iang))
                end if






            end do ! accumulate modes

            !
            ! Store value for point k.
            !
            output(iang) = s

        end do  ! evaluate points



        !
        ! Normalization
        !
        !output = output / real(nmodes,8)

    end subroutine idft_angle
    !******************************************************************












    !>  Computes the inverse discrete Fourier transform of an input data at 
    !!  user specified angles.
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
    subroutine idft_mode_angle(inreal, inimag, imode, angles, output)
        type(AD_D),       intent(in)      :: inreal
        type(AD_D),       intent(in)      :: inimag
        integer(ik),      intent(in)      :: imode
        real(rk),         intent(in)      :: angles(:)
        type(AD_D),       intent(inout)   :: output(:)

        type(AD_D)      :: sreal, simag, s, mag, phase
        real(rk)        :: freq
        integer         :: iang, nangles, quadrant


        nangles  = size(angles)


        !
        ! Initialize derivative arrays
        !
!        output = angles
        sreal  = inreal
        simag  = inreal
        s      = inreal
!        freq   = inreal

        output = ZERO
        sreal  = ZERO
        simag  = ZERO
        s      = ZERO
        freq   = ZERO




        !
        ! Loop over angles
        !
        do iang = 1,nangles

            !
            ! Accumulate contribution from each mode.
            !
            output(iang) = ZERO
            s            = ZERO

            
            !
            ! Compute current frequency
            !
            freq  = TWO * PI * real(imode-1,rk)

!            !
!            ! Compute magnitude and phase from the complex mode
!            !
!            print*, 'idft_mode_angle - inreal, inimag'
!            print*, inreal%xp_ad_
!            print*, inimag%xp_ad_
!            read*,
!            mag = sqrt( inreal**TWO + inimag**TWO )
!
!
!            if ( mag == ZERO ) then
!                phase = mag*ZERO ! avoid division by zero
!            else
!            
!                arg = inimag / mag
!
!                if ( arg == 1 ) then
!                    
!                phase = asin( inimag / mag )
!            end if
!
!
!
!            !
!            ! Check quadrant
!            !
!            if      ( (inreal >= ZERO)  .and. (inimag >= ZERO) ) then
!                quadrant = 1
!            else if ( (inreal < ZERO)   .and. (inimag >= ZERO) ) then
!                quadrant = 2
!            else if ( (inreal < ZERO)   .and. (inimag < ZERO)  ) then
!                quadrant = 3
!            else if ( (inreal >= ZERO)  .and. (inimag < ZERO)  ) then
!                quadrant = 4
!            end if
!
!
!            !
!            ! Account for quadrant
!            !
!            if      ( quadrant == 1 ) then
!                phase = phase
!            else if ( quadrant == 2) then
!                phase = PI - phase
!            else if ( quadrant == 3) then
!                phase = PI - phase
!            else if ( quadrant == 4) then
!                phase = phase
!            end if
!
!
!
!            print*, 'idft_mode_angle - s before'
!            print*, s%xp_ad_
!            print*, mag%xp_ad_
!            print*, phase%xp_ad_
!
!
!            if ( imode > 1 ) then
!                s = s +  TWO * mag*cos( freq * angles(iang) + phase )
!            else
!                s = s +  mag*cos( freq * angles(iang) + phase ) ! DC offset
!            end if
!
!            print*, 'idft_mode_angle - s after'
!            print*, s%xp_ad_
!            read*,



            if ( imode > 1 ) then
                s = s + TWO * (inreal*cos(freq*angles(iang)) - inimag*sin(freq*angles(iang)))
            else
                s = s + inreal*cos(freq*angles(iang)) - inimag*sin(freq*angles(iang))
            end if



            !
            ! Store value for point k.
            !
            output(iang) = s

        end do  ! evaluate points



        !
        ! Normalization
        !
        !output = output / real(nmodes,8)

    end subroutine idft_mode_angle
    !******************************************************************













    !>  Computes the inverse discrete Fourier transform of an input data at 
    !!  user specified points along the period, instead of the original evenly 
    !!  spaced periodic points.
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
    subroutine idft_points(ymin, periodicity, inreal, inimag, points, output)
        real(rk),         intent(in)  :: ymin
        real(rk),         intent(in)  :: periodicity
        type(AD_D),       intent(in)  :: inreal(:)
        type(AD_D),       intent(in)  :: inimag(:)
        type(AD_D),       intent(in)  :: points(:)
        type(AD_D),       intent(inout), allocatable  :: output(:)

        type(AD_D)   :: location(size(points))
        integer      :: npts, k

        npts    = size(points)

        !
        ! Compute location of each point in the period: 0-1
        !
        do k = 1,npts
            location(k) = (points(k)-ymin)/periodicity     ! 0-1.   0 = beginning of period.    1 = end of period.
        end do

        !
        ! Evaluate modes at specified locations
        !
        call idft_angle(inreal,inimag,location,output)

    end subroutine idft_points
    !******************************************************************












    !>  Computes the inverse discrete Fourier transform of an input data at 
    !!  user specified points along the period, instead of the original evenly 
    !!  spaced periodic points.
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
    subroutine idft_mode_points(ymin, periodicity, inreal, inimag, imode, points, output)
        real(rk),         intent(in)    :: ymin
        real(rk),         intent(in)    :: periodicity
        type(AD_D),       intent(in)    :: inreal
        type(AD_D),       intent(in)    :: inimag
        integer(ik),      intent(in)    :: imode
        real(rk),         intent(in)    :: points(:)
        type(AD_D),       intent(inout) :: output(:)

        real(rk)   :: angles(size(points))
        integer      :: npts, k

        npts = size(points)

        !
        ! Compute location of each point in the period: 0-1
        !
        do k = 1,npts
            angles(k) = (points(k)-ymin)/periodicity     ! 0-1.   0 = beginning of period.    1 = end of period.
        end do

        !
        ! Evaluate modes at specified locations
        !
        call idft_mode_angle(inreal,inimag,imode,angles,output)

    end subroutine idft_mode_points
    !******************************************************************





end module mod_dft
