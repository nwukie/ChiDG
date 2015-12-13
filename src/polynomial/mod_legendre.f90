module mod_legendre
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: XI_DIR,ETA_DIR,ZETA_DIR, &
                              ZERO, ONE, TWO, THREE, FOUR, FIVE, EIGHTH, HALF

    implicit none

contains


    !------------------------------------------------------------------
    !       polynomial values
    !------------------------------------------------------------------

    function LegendreVal(space_dim,currentmode,xpos,ypos,zpos) result(polyval)
        integer(kind=ik), intent(in)        :: space_dim,currentmode
        real(kind=rk),    intent(in)        :: xpos,ypos,zpos
        real(kind=rk)                       :: polyval

        select case (space_dim)
            case (1)    ! 1D
                polyval = LegendreVal1D(currentmode,xpos)
            case (2)    ! 2D
                polyval = LegendreVal2D(currentmode,xpos,ypos)
            case (3)    ! 3D
                polyval = LegendreVal3D(currentmode,xpos,ypos,zpos)
            case default
                print*, "Error - LegendreVal: valid space dimensions are (1,2,3)"
                stop
        end select
    end function




    recursive function LegendreVal1D(nterm,pos) result(polyval)
    !   Compute the value of the nterm
    !   legendre polynomial at the location pos
    !   between -1 and 1
    !   Edit list:  Nathan A. Wukie - 2/11/2015
        integer(ik),    intent(in) :: nterm
        real(rk),       intent(in) :: pos
        real(rk)                   :: polyval, polyval_nm1, polyval_nm2

        select case (nterm)
            ! Start recursion terms
            case (1)
                polyval = ONE
            case (2)
                polyval = pos
            case (3 :)
                ! Recursive definition for norder >= 2
                polyval_nm1=LegendreVal1D(nterm-1,pos)
                polyval_nm2=LegendreVal1D(nterm-2,pos)
                polyval = ((TWO*real(nterm-1,rk)-ONE)*pos*polyval_nm1 - &
                          ((real(nterm-1,rk)-ONE))*polyval_nm2)/real(nterm-1,rk)
        end select

!
!        select case (nterm)
!            case (1)
!                polyval = ONE
!            case (2)
!                polyval = pos
!            case (3)
!                polyval = HALF*(THREE*pos*pos - ONE)
!            case (4)
!                polyval = HALF*(FIVE*pos*pos*pos - THREE*pos)
!            case (5)
!                polyval = EIGHTH*(35._rk*pos*pos*pos*pos - 30._rk*pos*pos + THREE)
!            case (6)
!                polyval = EIGHTH*(63._rk*pos*pos*pos*pos*pos - 70._rk*pos*pos*pos + 15._rk*pos)
!            case default
!                stop "Error: Legendre polynomial - 1D mode limit exceeded"
!         end select

    end function



    function LegendreVal2D(currentmode,xi,eta) result(polyval)
    !   A set of 2D, Lagrange polynomials is associated with the coordinates
    !   in 'nodes'. This function computes the value of the Lagrange polynomial
    !   associated with the 'currentnode' at the coordinate '(xpos,ypos)'.
    !   Edit list:  Nathan A. Wukie - 2/15/2015
        use mod_ordering,   only: xi_order_2d,eta_order_2d

        integer(ik), intent(in)   :: currentmode
        real(rk),    intent(in)   :: xi,eta

        real(rk)                  :: polyval,xi_polyval,eta_polyval
        integer(ik)               :: xi_mode,eta_mode

        ! compute current x/y-node indices of 1D tensor product polynomials
        ! Example: for a 2D lagrange polynomial L(x,y)=L(x)*L(y), compute the
        ! x and y indices
        xi_mode  = xi_order_2d( currentmode)
        eta_mode = eta_order_2d(currentmode)

        xi_polyval  = LegendreVal1D(xi_mode,  xi)
        eta_polyval = LegendreVal1D(eta_mode,eta)

        polyval = xi_polyval*eta_polyval
    end function

    function LegendreVal3D(currentmode,xi,eta,zeta) result(polyval)
    !   A set of 3D, Lagrange polynomials is associated with the coordinates
    !   in 'nodes'. This function computes the value of the Lagrange polynomial
    !   associated with the 'currentnode' at the coordinate '(xpos,ypos,zpos)'.
    !   Edit list:  Nathan A. Wukie - 2/15/2015
        use mod_ordering,   only: xi_order_3d,eta_order_3d,zeta_order_3d

        integer(ik), intent(in)   :: currentmode
        real(rk),    intent(in)   :: xi,eta,zeta

        real(rk)                  :: polyval,xi_polyval,eta_polyval,zeta_polyval
        integer(ik)               :: xi_mode,eta_mode,zeta_mode

        ! compute current x/y-node indices of 1D tensor product polynomials
        ! Example: for a 2D lagrange polynomial L(x,y)=L(x)*L(y), compute the
        ! x and y indices
        xi_mode   = xi_order_3d(  currentmode)
        eta_mode  = eta_order_3d( currentmode)
        zeta_mode = zeta_order_3d(currentmode)

        xi_polyval   = LegendreVal1D(xi_mode,  xi)
        eta_polyval  = LegendreVal1D(eta_mode, eta)
        zeta_polyval = LegendreVal1D(zeta_mode,zeta)

        polyval = xi_polyval*eta_polyval*zeta_polyval
    end function










    !------------------------------------------------------------------
    !       polynomial derivatives
    !------------------------------------------------------------------

    function dlegendreVal(space_dim,currentmode,xi,eta,zeta,dir) result(dpolyval)
        integer(ik),   intent(in)          :: space_dim,currentmode
        real(rk),      intent(in)          :: xi,eta,zeta
        integer(ik),   intent(in)          :: dir

        real(rk)                           :: dpolyval

        select case (space_dim)
            case (1)    ! 1D
                dpolyval = DLegendreVal1D(currentmode,xi)
            case (2)    ! 2D
                dpolyval = DLegendreVal2D(currentmode,xi,eta,dir)
            case (3)    ! 3D
                dpolyval = DLegendreVal3D(currentmode,xi,eta,zeta,dir)
            case default
                print*, "Error - DLegendreVal: Valid space dimensions are (1,2,3)"
                stop
        end select
    end function


    recursive function DLegendreVal1D(nterm,pos) result(dpolyval)
    !   Compute the first derivative of the nterm
    !   Legendre polynomial at the location pos
    !   between -1 and 1
    !   Edit list:  Nathan A. Wukie - 2/11/2015
        integer(ik), intent(in)    :: nterm
        real(rk),    intent(in)    :: pos
        real(rk)                   :: dpolyval

        select case (nterm)
            ! Trivial evaluations
            case (1)
                dpolyval = ZERO
            case (2)
                dpolyval = ONE
            case (3 :)
                ! Recursive definition
                dpolyval = real(nterm-1,rk)*LegendreVal1D(nterm-1,pos) + pos*DLegendreVal1D(nterm-1,pos)

        end select
    end function






        function DLegendreVal2D(currentmode,xi,eta,dir) result(dpolyval)
    !   A set of 1D-Lagrange polynomials is associated with the coordinates
    !   in 'nodes'. This function compute the derivative of the Lagrange polynomial
    !   associated with the 'currentnode' at the location 'pos'.
    !   Edit list:  Nathan A. Wukie - 2/21/2015
        use mod_ordering,   only: xi_order_2d,eta_order_2d

        integer(ik), intent(in)            :: currentmode
        real(rk),    intent(in)            :: xi,eta
        integer(ik), intent(in)            :: dir

        integer(ik)                        :: xi_mode,eta_mode
        real(rk)                           :: dpolyval
        real(rk)                           :: xi_val,eta_val,dxi_val,deta_val

        xi_mode  = xi_order_2d(currentmode)
        eta_mode = eta_order_2d(currentmode)

        select case (dir)
            case (XI_DIR)
                dxi_val = DLegendreVal1D(xi_mode,xi)
                eta_val = LegendreVal1D(eta_mode,eta)

                dpolyval = dxi_val*eta_val
            case (ETA_DIR)

                xi_val   = LegendreVal1D(xi_mode,xi)
                deta_val = DLegendreVal1D(eta_mode,eta)

                dpolyval = xi_val*deta_val
            case default
                print*, "valid derivative directions are - 'xi', 'eta'"
                stop
        end select
    end function


    function DLegendreVal3D(currentmode,xi,eta,zeta,dir) result(dpolyval)
    !   A set of 3D-Lagrange polynomials is associated with the coordinates
    !   in 'nodes'. This function computes the derivative of the Lagrange polynomial
    !   associated with the 'currentnode' at the location 'pos' my multiplying the
    !   1D derivatives, since the 3D modes are constructed from a tensor product of 1D modes
    !   Edit list:  Nathan A. Wukie - 4/13/2015
        use mod_ordering,   only: xi_order_3d,eta_order_3d,zeta_order_3d

        integer(ik), intent(in)            :: currentmode
        real(rk),    intent(in)            :: xi,eta,zeta
        integer(ik), intent(in)            :: dir

        integer(ik)                        :: xi_mode,eta_mode,zeta_mode
        real(rk)                           :: dpolyval
        real(rk)                           :: xi_val,eta_val,zeta_val,dxi_val,deta_val,dzeta_val



        xi_mode   = xi_order_3d(  currentmode)
        eta_mode  = eta_order_3d( currentmode)
        zeta_mode = zeta_order_3d(currentmode)

        select case (dir)
            case (XI_DIR)
                dxi_val   = DLegendreVal1D(xi_mode,xi)
                eta_val   = LegendreVal1D(eta_mode,eta)
                zeta_val  = LegendreVal1D(zeta_mode,zeta)

                dpolyval  = dxi_val*eta_val*zeta_val
            case (ETA_DIR)
                xi_val    = LegendreVal1D(xi_mode,xi)
                deta_val  = DLegendreVal1D(eta_mode,eta)
                zeta_val  = LegendreVal1D(zeta_mode,zeta)

                dpolyval  = xi_val*deta_val*zeta_val
            case (ZETA_DIR)
                xi_val    = LegendreVal1D(xi_mode,xi)
                eta_val   = LegendreVal1D(eta_mode,eta)
                dzeta_val = DLegendreVal1D(zeta_mode,zeta)

                dpolyval  = xi_val*eta_val*dzeta_val
            case default
                print*, "valid derivative directions are - 'xi', 'eta', 'zeta'"
                stop
        end select
    end function









end module mod_legendre
