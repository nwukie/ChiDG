module mod_polynomial
#include <messenger.h>
    use mod_kinds,  only: rk,ik
    use type_point, only: point_t

    use mod_legendre
    use mod_lagrange
    use mod_ordering
    implicit none

contains




    !> Return value of a term in a polynomial expansion at a specified (xi, eta, zeta) location.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  spacedim    Integer specifying whether to compute the 1D, 2D, 
    !!                          or 3D expansion function. 
    !!  @param[in]  nterms      Number of terms in the polynomial expansion.
    !!  @param[in]  mode        Integer of the mode in the expansion to be computed.
    !!  @param[in]  node        point_t containing (xi,eta,zeta) location on reference element.
    !!
    !-----------------------------------------------------------------------------------------
    function PolynomialVal(spacedim,nterms,mode,node) result(polyval)
        integer(ik),    intent(in)  :: spacedim
        integer(ik),    intent(in)  :: nterms
        integer(ik),    intent(in)  :: mode
        type(point_t),  intent(in)  :: node

        real(rk)    :: xi, eta, zeta
        real(rk)    :: polyval


        xi   = node%c1_
        eta  = node%c2_
        zeta = node%c3_

        polyval = LegendreVal(spacedim,mode,xi,eta,zeta)

    end function PolynomialVal
    !*****************************************************************************************







    !>  Return directional derivative of a term in a polynomial expansion at a specified 
    !!  (xi, eta, zeta) location.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  spacedim   Integer specifying whether to compute the 1D, 2D, 
    !!                          or 3D expansion function.
    !!  @param[in]  nterms      Number of terms in the polynomial expansion.
    !!  @param[in]  mode        Integer of the mode in the expansion to be computed.
    !!  @param[in]  node        point_t containing (xi,eta,zeta) location on reference element.
    !!  @param[in]  dir         Integer indicating direction of the derivative (xi, eta, zeta).
    !!
    !----------------------------------------------------------------------------------------
    function DPolynomialVal(spacedim,nterms,mode,node,dir) result(dpolyval)
        integer(ik),    intent(in)  :: spacedim
        integer(ik),    intent(in)  :: nterms
        integer(ik),    intent(in)  :: mode
        type(point_t),  intent(in)  :: node
        integer(ik),    intent(in)  :: dir

        real(rk)    :: xi, eta, zeta
        real(rk)    :: dpolyval

        xi   = node%c1_
        eta  = node%c2_
        zeta = node%c3_

        dpolyval = DLegendreVal(spacedim,mode,xi,eta,zeta,dir)

    end function DPolynomialVal
    !****************************************************************************************





end module mod_polynomial
