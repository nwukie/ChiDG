module mod_cylindricalduct
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: HALF, ZERO, ONE, TWO
    use type_point,         only: point_t
    use type_function,      only: function_t
    use mod_function,       only: create_function
    use mod_gridspace,      only: linspace
    use mod_rootfinding,    only: bisect
    implicit none



contains



    !>  Compute eigenvalues of cylindrical duct modes for a given azimuthal
    !!  mode order 'm'.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/14/2016
    !!
    !!
    !---------------------------------------------------------------------------
    function compute_cylindricalduct_eigenvalues(m,nevalues) result(evalues)
        integer(ik),    intent(in)  :: m
        integer(ik),    intent(in)  :: nevalues

        real(rk), dimension(nevalues)   :: evalues
        real(rk), allocatable           :: xvals(:)
        real(rk)                        :: fnew, fold, xnew, xold
        type(point_t)                   :: x
        integer(ik)                     :: ieig, ix
        class(function_t), allocatable  :: efcn


        !
        ! Create eigenfunction and set 'm'
        !
        call create_function(efcn,'cylindricalduct_eigenfunction')
        call efcn%set_option('m',real(m,rk))

        !
        ! Initialize values
        !
        evalues = ZERO
        call x%set(ZERO,ZERO,ZERO)


        !
        ! Numerically compute each eigen value
        !
        do ieig = 1,nevalues

            !
            ! Eigenvalue for (0,1) is ZERO. Bisection
            ! wouldn't find this because iteration starts 
            ! after 0.
            !
            if ( m == 0 ) then
                if ( ieig == 1 ) then
                    evalues(ieig) = ZERO
                    exit
                end if
            end if


            
            !
            ! Get discretization values to try
            !
            if (ieig == 1) then
                ! Start from the beginning of the function
                xvals = linspace(0.1_rk,100._rk,1000)
            else
                ! Start after the previous zero
                xvals = linspace(evalues(ieig-1)+0.1,100._rk,10000)
            end if


            ! 
            ! Iterate through xvals and look for sign change
            !
            do ix = 1,size(xvals)
                x%c1_ = xvals(ix)
                fnew = efcn%compute(ZERO,x)
                xnew = xvals(ix)

                if (ix > 1) then
                    ! Test for sign change
                    if ( int(sign(ONE,fnew)) /= int(sign(ONE,fold)) ) then
                        ! Sign changed, so start bisection to find root.
                        evalues(ieig) = bisect(efcn, xold, xnew)
                        exit
                    else
                        fold = fnew
                        xold = xnew
                    end if
                else
                    fold = fnew
                    xold = xnew
                end if

            end do




        end do





    end function compute_cylindricalduct_eigenvalues
    !***************************************************************************




















    !>  Return the value of a specified cylindrical duct mode.
    !!
    !!  For a specified azimuthal mode(m), and the eigenvalue of an associated 
    !!  radial mode(alpha), compute the value of the mode at radial location(r),
    !!  bounded by the outer wall(rmax)
    !!
    !!  The mode values are normalized by N such that:
    !!
    !!  \f$ \int_0^1 N_{m \mu}^2 J(\alpha_{m \mu})^2 r dr = 1 \f$
    !!
    !!  This is outlined in the document:
    !!  "Fundamentals of Duct Acoustics" by Sjoerd W. Rienstra
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/14/2016
    !!
    !----------------------------------------------------------------------------
    elemental function compute_cylindricalduct_mode(m,alpha,r,rmax) result(val)
        integer(ik),    intent(in)  :: m
        real(rk),       intent(in)  :: alpha
        real(rk),       intent(in)  :: r
        real(rk),       intent(in)  :: rmax

        real(rk)    :: val, normalization


        ! Compute value of mode at location 'r'
        val = bessel_jn(m,alpha*r/rmax)


        ! Compute normalization
        if ( m == 0 ) then
            normalization = sqrt(TWO)
        else
            normalization = sqrt(TWO)/(bessel_jn(m,alpha)*sqrt(ONE - (real(m,rk)**TWO)/alpha**TWO))
        end if

        ! Apply normalization 
        val = val*normalization

    end function compute_cylindricalduct_mode
    !*****************************************************************************




























end module mod_cylindricalduct
