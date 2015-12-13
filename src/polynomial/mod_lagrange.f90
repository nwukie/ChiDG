module mod_lagrange
    use mod_kinds,      only: rk,ik
    use mod_constants
    use mod_GaussLegendre, only: gl_nodes

    implicit none

contains

    !------------------------------------------------------------------
    !       polynomial values
    !------------------------------------------------------------------

    function LagrangeVal(space_dim,nterms,currentnode,xpos,ypos,zpos) result(polyval)
        integer(kind=ik), intent(in)        :: space_dim,nterms,currentnode
        real(kind=rk), intent(in)           :: xpos,ypos,zpos
        real(kind=rk)                       :: polyval

        select case (space_dim)
            case (1)    ! 1D
                polyval = LagrangeVal1D(nterms,currentnode,xpos)
            case (2)    ! 2D
                polyval = LagrangeVal2D(nterms,currentnode,xpos,ypos)
            case (3)    ! 3D
                polyval = LagrangeVal3D(nterms,currentnode,xpos,ypos,zpos)
            case default
                print*, "Error - lagrangeVal: valid space dimensions are (1,2,3)"
                stop
        end select
    end function

    function LagrangeVal1D(nterms,currentnode,pos) result(polyval)
    !   A set of 1D, Lagrange polynomials is associated with the coordinates
    !   in 'nodes'. This function computes the value of the Lagrange polynomial
    !   associated with the 'currentnode' at the coordinate 'pos'.
    !   Edit list:  Nathan A. Wukie - 2/15/2015
        integer(kind=ik), intent(in)                :: nterms,currentnode
        real(kind=rk), intent(in)                   :: pos
        real(kind=rk), dimension(nterms)            :: nodes
        real(kind=rk)                               :: polyval
        integer(kind=ik)                            :: j
        real(kind=rk)                               :: xi,xj

        ! Get Gauss-Legendre nodes
        call gl_nodes(nterms,nodes)

        polyval = 1._rk
        xi = nodes(currentnode)
        select case (nterms)
            case (1)
                polyval = 1._rk
            case (2 :)
                do j=1,nterms
                    xj = nodes(j)
                    if (currentnode /= j) then
                        polyval = polyval*(pos-xj)/(xi-xj)
                    end if
                end do
        end select

    end function

    function LagrangeVal2D(nterms,currentmode,xi,eta) result(polyval)
    !   A set of 2D, Lagrange polynomials is associated with the coordinates
    !   in 'nodes'. This function computes the value of the Lagrange polynomial
    !   associated with the 'currentnode' at the coordinate '(xpos,ypos)'.
    !   Edit list:  Nathan A. Wukie - 2/15/2015
        use mod_ordering,   only: xi_order_2d,eta_order_2d

        integer(kind=ik), intent(in)                :: nterms,currentmode
        real(kind=rk),    intent(in)                :: xi,eta

        real(kind=rk)                               :: polyval,xi_polyval,eta_polyval
        integer(kind=ik)                            :: nterms1d,xi_mode,eta_mode

        ! compute the number of terms in the 1D polynomial
        ! really just finding the square root of nterms, but there isn't
        ! such a function for integers in fortran
        nterms1d = 0
        do while (nterms1d*nterms1d /= nterms)
            nterms1d = nterms1d+1
        end do

        ! compute current x/y-node indices of 1D tensor product polynomials
        ! Example: for a 2D lagrange polynomial L(x,y)=L(x)*L(y), compute the
        ! x and y indices
!        print*, currentmode
        xi_mode  = xi_order_2d( currentmode)
        eta_mode = eta_order_2d(currentmode)

        xi_polyval  = LagrangeVal1D(nterms1d,xi_mode,  xi)
        eta_polyval = LagrangeVal1D(nterms1d,eta_mode,eta)

        polyval = xi_polyval*eta_polyval
    end function



    function LagrangeVal3D(nterms,currentmode,xi,eta,zeta) result(polyval)
    !   A set of 3D, Lagrange polynomials is associated with the coordinates
    !   in 'nodes'. This function computes the value of the Lagrange polynomial
    !   associated with the 'currentnode' at the coordinate '(xpos,ypos,zpos)'.
    !   Edit list:  Nathan A. Wukie - 2/15/2015
        use mod_ordering,   only: xi_order_3d,eta_order_3d,zeta_order_3d

        integer(kind=ik), intent(in)                :: nterms,currentmode
        real(kind=rk),    intent(in)                :: xi,eta,zeta

        real(kind=rk)                               :: polyval,xi_polyval,eta_polyval,zeta_polyval
        integer(kind=ik)                            :: nterms1d,xi_mode,eta_mode,zeta_mode

        ! compute the number of terms in the 1D polynomial
        ! really just finding the cubed-root of nterms, but there isn't
        ! such a function for integers in fortran
        nterms1d = 0
        do while (nterms1d*nterms1d*nterms1d /= nterms)
            nterms1d = nterms1d+1
        end do

        ! compute current x/y-node indices of 1D tensor product polynomials
        ! Example: for a 2D lagrange polynomial L(x,y)=L(x)*L(y), compute the
        ! x and y indices
        xi_mode   = xi_order_3d(  currentmode)
        eta_mode  = eta_order_3d( currentmode)
        zeta_mode = zeta_order_3d(currentmode)

        xi_polyval   = LagrangeVal1D(nterms1d,xi_mode,  xi)
        eta_polyval  = LagrangeVal1D(nterms1d,eta_mode, eta)
        zeta_polyval = LagrangeVal1D(nterms1d,zeta_mode,zeta)

        polyval = xi_polyval*eta_polyval*zeta_polyval
    end function


    !------------------------------------------------------------------
    !       polynomial derivatives
    !------------------------------------------------------------------

    function DLagrangeVal(space_dim,nterms,currentmode,xi,eta,zeta,dir) result(dpolyval)
        integer(kind=ik),   intent(in)          :: space_dim,nterms,currentmode
        real(kind=rk),      intent(in)          :: xi,eta,zeta
        integer(kind=ik),   intent(in)          :: dir
!        character(len=*),   intent(in)          :: dir

        real(kind=rk)                           :: dpolyval

        select case (space_dim)
            case (1)    ! 1D
                dpolyval = DLagrangeVal1D(nterms,currentmode,xi)
            case (2)    ! 2D
                dpolyval = DLagrangeVal2D(nterms,currentmode,xi,eta,dir)
            case (3)    ! 3D
                dpolyval = DLagrangeVal3D(nterms,currentmode,xi,eta,zeta,dir)
            case default
                print*, "Error - DLagrangeVal: Valid space dimensions are (1,2,3)"
                stop
        end select
    end function


    function DLagrangeVal1D(nterms1d,currentmode,xi) result(dpolyval)
    !   A set of 1D-Lagrange polynomials is associated with the coordinates
    !   in 'nodes'. This function compute the derivative of the Lagrange polynomial
    !   associated with the 'currentnode' at the location 'pos'.
    !   Edit list:  Nathan A. Wukie - 2/21/2015

        integer(kind=ik), intent(in)            :: nterms1d,currentmode
        real(kind=rk),    intent(in)            :: xi

        real(kind=rk), dimension(nterms1d)      :: nodes
        integer(kind=ik)                        :: j,k
        real(kind=rk)                           :: dpolyval
        real(kind=rk)                           :: prod,xi_i,xi_j,xk

        call gl_nodes(nterms1d,nodes)

        dpolyval = 0._rk
        xi_i = nodes(currentmode)
        do j=1,nterms1d
            xi_j = nodes(j)
            prod = 1._rk
            if (j /= currentmode) then
                do k=1,nterms1d
                    xk = nodes(k)
                    if ( k /= currentmode .and. k /= j) then
                        prod = prod*(xi-xk)/(xi_i-xk)
                    end if
                end do
                dpolyval = dpolyval + prod/(xi_i-xi_j)
            end if
        end do

    end function

    function DLagrangeVal2D(nterms,currentmode,xi,eta,dir) result(dpolyval)
    !   A set of 1D-Lagrange polynomials is associated with the coordinates
    !   in 'nodes'. This function compute the derivative of the Lagrange polynomial
    !   associated with the 'currentnode' at the location 'pos'.
    !   Edit list:  Nathan A. Wukie - 2/21/2015
        use mod_ordering,   only: xi_order_2d,eta_order_2d

        integer(kind=ik), intent(in)            :: nterms,currentmode
        real(kind=rk),    intent(in)            :: xi,eta
        integer(kind=ik), intent(in)            :: dir
!        character(len=*), intent(in)            :: dir

        integer(kind=ik)                        :: nterms1d,xi_mode,eta_mode
        real(kind=rk)                           :: dpolyval
        real(kind=rk)                           :: xi_val,eta_val,dxi_val,deta_val

        ! compute the number of terms in the 1d polynomial
        ! really just finding the square root of nterms, but there isn't
        ! such a function for integers in fortran
        nterms1d = 0
        do while (nterms1d*nterms1d /= nterms)
            nterms1d = nterms1d+1
        end do

        xi_mode  = xi_order_2d(currentmode)
        eta_mode = eta_order_2d(currentmode)

        select case (dir)
            case (XI_DIR)
                dxi_val = DLagrangeVal1D(nterms1d,xi_mode,xi)
                eta_val = LagrangeVal1D(nterms1d,eta_mode,eta)

                dpolyval = dxi_val*eta_val
            case (ETA_DIR)

                xi_val   = LagrangeVal1D(nterms1d,xi_mode,xi)
                deta_val = DLagrangeVal1D(nterms1d,eta_mode,eta)

                dpolyval = xi_val*deta_val
            case default
                write(*,*) "valid derivative directions are - 'xi', 'eta'"
                stop
        end select
    end function


    function DLagrangeVal3D(nterms,currentmode,xi,eta,zeta,dir) result(dpolyval)
    !   A set of 3D-Lagrange polynomials is associated with the coordinates
    !   in 'nodes'. This function computes the derivative of the Lagrange polynomial
    !   associated with the 'currentnode' at the location 'pos' my multiplying the
    !   1D derivatives, since the 3D modes are constructed from a tensor product of 1D modes
    !   Edit list:  Nathan A. Wukie - 4/13/2015
        use mod_ordering,   only: xi_order_3d,eta_order_3d,zeta_order_3d

        integer(kind=ik), intent(in)            :: nterms,currentmode
        real(kind=rk),    intent(in)            :: xi,eta,zeta
        integer(kind=ik), intent(in)            :: dir

        integer(kind=ik)                        :: nterms1d,xi_mode,eta_mode,zeta_mode
        real(kind=rk)                           :: dpolyval
        real(kind=rk)                           :: xi_val,eta_val,zeta_val,dxi_val,deta_val,dzeta_val

        ! compute the number of terms in the 1d polynomial
        ! really just finding the square root of nterms, but there isn't
        ! such a function for integers in fortran
        nterms1d = 0
        do while (nterms1d*nterms1d*nterms1d /= nterms)
            nterms1d = nterms1d+1
        end do

        xi_mode   = xi_order_3d(  currentmode)
        eta_mode  = eta_order_3d( currentmode)
        zeta_mode = zeta_order_3d(currentmode)

        select case (dir)
            case (XI_DIR)
                dxi_val   = DLagrangeVal1D(nterms1d,xi_mode,xi)
                eta_val   = LagrangeVal1D(nterms1d,eta_mode,eta)
                zeta_val  = LagrangeVal1D(nterms1d,zeta_mode,zeta)

                dpolyval  = dxi_val*eta_val*zeta_val
            case (ETA_DIR)
                xi_val    = LagrangeVal1D(nterms1d,xi_mode,xi)
                deta_val  = DLagrangeVal1D(nterms1d,eta_mode,eta)
                zeta_val  = LagrangeVal1D(nterms1d,zeta_mode,zeta)

                dpolyval  = xi_val*deta_val*zeta_val
            case (ZETA_DIR)
                xi_val    = LagrangeVal1D(nterms1d,xi_mode,xi)
                eta_val   = LagrangeVal1D(nterms1d,eta_mode,eta)
                dzeta_val = DLagrangeVal1D(nterms1d,zeta_mode,zeta)

                dpolyval  = xi_val*eta_val*dzeta_val
            case default
                write(*,*) "valid derivative directions are - 'xi', 'eta', 'zeta'"
                stop
        end select
    end function



end module mod_lagrange
