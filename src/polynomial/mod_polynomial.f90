module mod_polynomial
    use mod_kinds,  only: rk,ik

    use type_point, only: point_t

    use mod_legendre
    use mod_lagrange

    implicit none

contains

    !>
    !!
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------
    function PolynomialVal(space_dim,nterms,mode,node) result(polyval)
        integer(ik),    intent(in)   :: space_dim,nterms,mode
        type(point_t),  intent(in)   :: node

        real(rk)                     :: xi,eta,zeta
        real(rk)                     :: polyval

        xi   = node%c1_
        eta  = node%c2_
        zeta = node%c3_

        polyval = LegendreVal(space_dim,mode,xi,eta,zeta)

    end function
    !###############################################################################




    !>
    !!
    !!
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------
    function DPolynomialVal(space_dim,nterms,mode,node,dir) result(dpolyval)
        integer(ik),    intent(in)   :: space_dim,nterms,mode
        type(point_t),  intent(in)   :: node
        integer(ik),    intent(in)   :: dir

        real(rk)                     :: xi,eta,zeta
        real(rk)                     :: dpolyval

        xi   = node%c1_
        eta  = node%c2_
        zeta = node%c3_

        dpolyval = DLegendreVal(space_dim,mode,xi,eta,zeta,dir)

    end function
    !###############################################################################





end module mod_polynomial
