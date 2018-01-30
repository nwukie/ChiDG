module mod_legendre
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: XI_DIR,ETA_DIR,ZETA_DIR, &
                              ZERO, ONE, TWO, THREE, FOUR, FIVE, EIGHTH, HALF
    use mod_ordering,   only: xi_order_2d, eta_order_2d, &
                              xi_order_3d, eta_order_3d, zeta_order_3d
    use ieee_arithmetic,    only: ieee_is_nan

    implicit none

contains


    !> Compute value of hierarchical Legendre polynomial expansion.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/20/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function legendre_val(space_dim,currentmode,xpos,ypos,zpos) result(polyval)
        integer(ik),    intent(in)  :: space_dim
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xpos
        real(rk),       intent(in)  :: ypos
        real(rk),       intent(in)  :: zpos

        real(rk)                    :: polyval

        select case (space_dim)
            case (1)    ! 1D
                polyval = legendre_val1D(currentmode,xpos)
            case (2)    ! 2D
                polyval = legendre_val2D(currentmode,xpos,ypos)
            case (3)    ! 3D
                polyval = legendre_val3D(currentmode,xpos,ypos,zpos)
            case default
                print*, "Error - legendre_val: valid space dimensions are (1,2,3)"
                stop
        end select

    end function legendre_val
    !*****************************************************************************************




    



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/20/2016
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    recursive function legendre_val1D(nterm,pos) result(polyval)
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
                polyval_nm1=legendre_val1D(nterm-1,pos)
                polyval_nm2=legendre_val1D(nterm-2,pos)
                polyval = ((TWO*real(nterm-1,rk)-ONE)*pos*polyval_nm1 - &
                          ((real(nterm-1,rk)-ONE))*polyval_nm2)/real(nterm-1,rk)
        end select


    end function legendre_val1D
    !****************************************************************************************











    !>  A set of 2D, Lagrange polynomials is associated with the coordinates
    !!  in 'nodes'. This function computes the value of the Lagrange polynomial
    !!  associated with the 'currentnode' at the coordinate '(xpos,ypos)'.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   3/20/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function legendre_val2D(currentmode,xi,eta) result(polyval)
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xi
        real(rk),       intent(in)  ::eta

        real(rk)    :: polyval, xi_polyval, eta_polyval
        integer(ik) :: xi_mode, eta_mode

        ! compute current x/y-node indices of 1D tensor product polynomials
        ! Example: for a 2D lagrange polynomial L(x,y)=L(x)*L(y), compute the
        ! x and y indices
        xi_mode  = xi_order_2d( currentmode)
        eta_mode = eta_order_2d(currentmode)

        xi_polyval  = legendre_val1D(xi_mode,  xi)
        eta_polyval = legendre_val1D(eta_mode,eta)

        polyval = xi_polyval*eta_polyval

    end function legendre_val2D
    !*****************************************************************************************










    !>  A set of 3D, Lagrange polynomials is associated with the coordinates
    !!  in 'nodes'. This function computes the value of the Lagrange polynomial
    !!  associated with the 'currentnode' at the coordinate '(xpos,ypos,zpos)'.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   3/20/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function legendre_val3D(currentmode,xi,eta,zeta) result(polyval)
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xi
        real(rk),       intent(in)  :: eta
        real(rk),       intent(in)  :: zeta

        real(rk)    :: polyval, xi_polyval, eta_polyval, zeta_polyval
        integer(ik) :: xi_mode, eta_mode, zeta_mode

        ! compute current x/y-node indices of 1D tensor product polynomials
        ! Example: for a 2D lagrange polynomial L(x,y)=L(x)*L(y), compute the
        ! x and y indices
        xi_mode   = xi_order_3d(  currentmode)
        eta_mode  = eta_order_3d( currentmode)
        zeta_mode = zeta_order_3d(currentmode)

        xi_polyval   = legendre_val1D(xi_mode,  xi)
        eta_polyval  = legendre_val1D(eta_mode, eta)
        zeta_polyval = legendre_val1D(zeta_mode,zeta)

        polyval = xi_polyval*eta_polyval*zeta_polyval

    end function legendre_val3D
    !*****************************************************************************************










    !> Compute directional derivative of legendre polynomial.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/20/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function dlegendre_val(space_dim,currentmode,xi,eta,zeta,dir) result(dpolyval)
        integer(ik),    intent(in)  :: space_dim
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xi
        real(rk),       intent(in)  :: eta
        real(rk),       intent(in)  :: zeta
        integer(ik),    intent(in)  :: dir

        real(rk)    :: dpolyval

        select case (space_dim)
            case (1)    ! 1D
                dpolyval = dlegendre_val1D(currentmode,xi)
            case (2)    ! 2D
                dpolyval = dlegendre_val2D(currentmode,xi,eta,dir)
            case (3)    ! 3D
                dpolyval = dlegendre_val3D(currentmode,xi,eta,zeta,dir)
            case default
                print*, "Error - dlegendre_val: Valid space dimensions are (1,2,3)"
                stop
        end select

    end function dlegendre_val
    !*****************************************************************************************







    !>  Compute the first derivative of the nterm Legendre polynomial at the location 
    !!  'pos' between -1 and 1.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/20/2016
    !!
    !-----------------------------------------------------------------------------------------
    recursive function dlegendre_val1D(nterm,pos) result(dpolyval)
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
                dpolyval = real(nterm-1,rk)*legendre_val1D(nterm-1,pos) + pos*dlegendre_val1D(nterm-1,pos)

        end select

    end function dlegendre_val1D
    !*****************************************************************************************












    !>  A set of 1D-Lagrange polynomials is associated with the coordinates
    !!  in 'nodes'. This function compute the derivative of the Lagrange polynomial
    !!  associated with the 'currentnode' at the location 'pos'.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/20/2016
    !!
    !-----------------------------------------------------------------------------------------
    function dlegendre_val2D(currentmode,xi,eta,dir) result(dpolyval)
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xi
        real(rk),       intent(in)  :: eta
        integer(ik),    intent(in)  :: dir

        integer(ik) :: xi_mode, eta_mode
        real(rk)    :: dpolyval
        real(rk)    :: xi_val, eta_val, dxi_val, deta_val

        xi_mode  = xi_order_2d(currentmode)
        eta_mode = eta_order_2d(currentmode)

        select case (dir)
            case (XI_DIR)
                dxi_val = dlegendre_val1D(xi_mode,xi)
                eta_val = legendre_val1D(eta_mode,eta)

                dpolyval = dxi_val*eta_val

            case (ETA_DIR)
                xi_val   = legendre_val1D(xi_mode,xi)
                deta_val = dlegendre_val1D(eta_mode,eta)

                dpolyval = xi_val*deta_val

            case (ZETA_DIR)
                ! By definition of 2D polynomial, no derivative in ZETA dimension
                dpolyval = ZERO 

            case default
                print*, "valid derivative directions are - 'XI_DIR', 'ETA_DIR', 'ZETA_DIR'"
                stop
        end select

    end function dlegendre_val2D
    !*****************************************************************************************






    !>  A set of 3D-Lagrange polynomials is associated with the coordinates
    !!  in 'nodes'. This function computes the derivative of the Lagrange polynomial
    !!  associated with the 'currentnode' at the location 'pos' my multiplying the
    !!  1D derivatives, since the 3D modes are constructed from a tensor product of 1D modes
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/20/2016
    !!
    !----------------------------------------------------------------------------------------
    function dlegendre_val3D(currentmode,xi,eta,zeta,dir) result(dpolyval)
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xi
        real(rk),       intent(in)  :: eta
        real(rk),       intent(in)  :: zeta
        integer(ik),    intent(in)  :: dir

        integer(ik) :: xi_mode, eta_mode, zeta_mode
        real(rk)    :: dpolyval
        real(rk)    :: xi_val, eta_val, zeta_val, dxi_val, deta_val, dzeta_val



        xi_mode   = xi_order_3d(  currentmode)
        eta_mode  = eta_order_3d( currentmode)
        zeta_mode = zeta_order_3d(currentmode)

        select case (dir)
            case (XI_DIR)
                dxi_val   = dlegendre_val1D(xi_mode,xi)
                eta_val   = legendre_val1D(eta_mode,eta)
                zeta_val  = legendre_val1D(zeta_mode,zeta)

                dpolyval  = dxi_val*eta_val*zeta_val
            case (ETA_DIR)
                xi_val    = legendre_val1D(xi_mode,xi)
                deta_val  = dlegendre_val1D(eta_mode,eta)
                zeta_val  = legendre_val1D(zeta_mode,zeta)

                dpolyval  = xi_val*deta_val*zeta_val
            case (ZETA_DIR)
                xi_val    = legendre_val1D(xi_mode,xi)
                eta_val   = legendre_val1D(eta_mode,eta)
                dzeta_val = dlegendre_val1D(zeta_mode,zeta)

                dpolyval  = xi_val*eta_val*dzeta_val
            case default
                print*, "valid derivative directions are - 'XI_DIR', 'ETA_DIR', 'ZETA_DIR'"
                stop
        end select

    end function dlegendre_val3D
    !*****************************************************************************************










    !>  Second/mixed derivatives of modes in legendre tensor product basis.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   8/7/2017
    !!
    !-----------------------------------------------------------------------------------------
    function ddlegendre_val(space_dim,currentmode,xi,eta,zeta,partial1,partial2) result(res)
        integer(ik),    intent(in)  :: space_dim
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xi
        real(rk),       intent(in)  :: eta
        real(rk),       intent(in)  :: zeta
        integer(ik),    intent(in)  :: partial1
        integer(ik),    intent(in)  :: partial2

        real(rk)    :: res

        select case (space_dim)
            case (1)    ! 1D
                res = ddlegendre_val1D(currentmode,xi)
            case (2)    ! 2D
                res = ddlegendre_val2D(currentmode,xi,eta,partial1,partial2)
            case (3)    ! 3D
                res = ddlegendre_val3D(currentmode,xi,eta,zeta,partial1,partial2)
            case default
                print*, "Error - ddlegendre_val: Valid space dimensions are (1,2,3)"
                stop
        end select

    end function ddlegendre_val
    !*****************************************************************************************














    !>  Compute the second derivative of the nterm Legendre polynomial at the location 
    !!  'pos' between -1 and 1.
    !!
    !!  The Legendre polynomials satisfy the Legendre equation:
    !!      ddx( (1-x^2)*ddx(Pn) )  +  n*(n+1)*Pn = 0
    !!
    !!  Applying the chain rule:
    !!      ddx(1 - x^2) * ddx(Pn)  -  2*x* ddx(ddx(Pn))  +  n*(n+1)*Pn = 0
    !!
    !!  Rearranging to solve for the second derivative: ddx(ddx(Pn))
    !!      ddx(ddx(Pn)) = [2*x*ddx(Pn)  -  n*(n+1)*Pn]/(1-x^2)
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   8/07/2017
    !!
    !-----------------------------------------------------------------------------------------
    recursive function ddlegendre_val1d(nterm,pos) result(res)
        use mod_constants,  only: ZERO, ONE, TWO, THREE, SEVEN, EIGHT
        integer(ik), intent(in)    :: nterm
        real(rk),    intent(in)    :: pos

        real(rk)    :: res

        ! DANGEROUS AT (pos=1)
        !res = (TWO*pos*dlegendre_val1D(nterm,pos) - real(nterm-1,rk)*(real(nterm-1,rk) + ONE)*legendre_val1D(nterm,pos))/(ONE-pos*pos)

        select case(nterm)
            case(1)
                res = ZERO
            case(2)
                res = ZERO
            case(3)
                res = THREE
            case(4) 
                res = 15._rk*pos
            case(5)
                res = (15._rk/TWO)*(SEVEN*pos*pos - ONE)
            case(6)
                res = (105_rk/TWO)*pos*(THREE*pos*pos - ONE)
            case(7)
                res = (105_rk/EIGHT)*(33._rk*pos*pos*pos*pos - 18._rk*pos*pos + ONE)
            case(8)
                res = (63._rk/EIGHT)*pos*(143._rk*pos*pos*pos*pos - 110._rk*pos*pos + 15._rk)
            case default
                print*, "Error: ddlegendre_val1d is only defined through order 8"
                stop
        end select

    end function ddlegendre_val1d
    !*****************************************************************************************




    !>  A set of 1D-Lagrange polynomials is associated with the coordinates
    !!  in 'nodes'. This function compute the derivative of the 2D Legendre modal tensor 
    !!  product associated with the 'currentnode' at the location 'pos'.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   8/7/2017
    !!
    !-----------------------------------------------------------------------------------------
    function ddlegendre_val2d(currentmode,xi,eta,partial1,partial2) result(res)
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xi
        real(rk),       intent(in)  :: eta
        integer(ik),    intent(in)  :: partial1
        integer(ik),    intent(in)  :: partial2

        integer(ik) :: xi_mode, eta_mode
        real(rk)    :: term1, term2, res

        !
        ! Check valid input for derivatives
        !
        if ((partial1 /= XI_DIR) .and. (partial1 /= ETA_DIR) .and. (partial1 /= ZETA_DIR) .or. &
            (partial2 /= XI_DIR) .and. (partial2 /= ETA_DIR) .and. (partial2 /= ZETA_DIR) ) then
            print*, "valid derivative directions are - 'XI_DIR', 'ETA_DIR', 'ZETA_DIR'"
            stop
        end if


        !
        ! Get indices of terms that construct tensor product mode
        !
        xi_mode  = xi_order_2d(currentmode)
        eta_mode = eta_order_2d(currentmode)


        !
        ! Pure second derivative
        !
        if (partial1 == partial2) then

            select case (partial1)
                case (XI_DIR)
                    term1 = ddlegendre_val1D(xi_mode,xi)
                    term2 =   legendre_val1D(eta_mode,eta)

                case (ETA_DIR)
                    term1 =   legendre_val1D(xi_mode,xi)
                    term2 = ddlegendre_val1D(eta_mode,eta)
            end select

        !
        ! Mixed derivative:  dd(phi)/dxideta = [d(phi)/dxi][d(phi)/deta]
        !
        else
            term1 = dlegendre_val1D(xi_mode,xi)
            term2 = dlegendre_val1D(eta_mode,eta)

        end if


        res = term1 * term2


        ! By definition of 2D polynomial, no derivative in ZETA dimension
        if ((partial1 == ZETA_DIR) .or. (partial2 == ZETA_DIR)) res = ZERO


    end function ddlegendre_val2d
    !*****************************************************************************************






    !>  A set of 1D-Legendre polynomials is associated with the coordinates
    !!  in 'nodes'. This function compute the derivative of the 3D Legendre modal tensor 
    !!  product associated with the 'currentnode' at the location 'pos'.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   8/7/2016
    !!
    !-----------------------------------------------------------------------------------------
    function ddlegendre_val3d(currentmode,xi,eta,zeta,partial1,partial2) result(res)
        integer(ik),    intent(in)  :: currentmode
        real(rk),       intent(in)  :: xi
        real(rk),       intent(in)  :: eta
        real(rk),       intent(in)  :: zeta
        integer(ik),    intent(in)  :: partial1
        integer(ik),    intent(in)  :: partial2

        integer(ik) :: xi_mode, eta_mode, zeta_mode
        real(rk)    :: term1, term2, term3, res

        xi_mode   = xi_order_3d(currentmode)
        eta_mode  = eta_order_3d(currentmode)
        zeta_mode = zeta_order_3d(currentmode)


        !
        ! Check valid input for derivatives
        !
        if ((partial1 /= XI_DIR) .and. (partial1 /= ETA_DIR) .and. (partial1 /= ZETA_DIR) .or. &
            (partial2 /= XI_DIR) .and. (partial2 /= ETA_DIR) .and. (partial2 /= ZETA_DIR) ) then
            print*, "valid derivative directions are - 'XI_DIR', 'ETA_DIR', 'ZETA_DIR'"
            stop
        end if




        !
        ! Pure second derivative
        !
        if (partial1 == partial2) then

            select case (partial1)
                case (XI_DIR)
                    term1 = ddlegendre_val1D(xi_mode,xi)
                    term2 =   legendre_val1D(eta_mode,eta)
                    term3 =   legendre_val1D(zeta_mode,zeta)

                case (ETA_DIR)
                    term1 =   legendre_val1D(xi_mode,xi)
                    term2 = ddlegendre_val1D(eta_mode,eta)
                    term3 =   legendre_val1D(zeta_mode,zeta)

                case (ZETA_DIR)
                    term1 =   legendre_val1D(xi_mode,xi)
                    term2 =   legendre_val1D(eta_mode,eta)
                    term3 = ddlegendre_val1D(zeta_mode,zeta)

            end select


        !
        ! Mixed derivative:  
        !   ex:  dd(phi)/dxideta   = [d(phi)/dxi] * [d(phi)/deta] * [phi(zeta)]
        !   ex:  dd(phi)/detadzeta = [phi(xi)]  *  [d(phi)/deta] * [d(phi)/dzeta]
        !
        else


            select case(partial1)
                case (XI_DIR)
                    term1 = dlegendre_val1D(xi_mode,xi)
                case (ETA_DIR)
                    term1 = dlegendre_val1D(eta_mode,eta)
                case (ZETA_DIR)
                    term1 = dlegendre_val1D(zeta_mode,zeta)
            end select


            select case(partial2)
                case (XI_DIR)
                    term2 = dlegendre_val1D(xi_mode,xi)
                case (ETA_DIR)
                    term2 = dlegendre_val1D(eta_mode,eta)
                case (ZETA_DIR)
                    term2 = dlegendre_val1D(zeta_mode,zeta)
            end select

            
            ! Determine third term by knowledge of first two terms in derivative
            select case(partial1 + partial2)
                case (XI_DIR + ETA_DIR)
                    term3 = legendre_val1D(zeta_mode,zeta)
                case (XI_DIR + ZETA_DIR)
                    term3 = legendre_val1D(eta_mode,eta)
                case (ETA_DIR + ZETA_DIR)
                    term3 = legendre_val1D(xi_mode,xi)
            end select



        end if


        ! Compute output
        res = term1 * term2 * term3

        if (ieee_is_nan(res)) print*, 'term ', currentmode, ' is nan'


    end function ddlegendre_val3d
    !*****************************************************************************************










end module mod_legendre
