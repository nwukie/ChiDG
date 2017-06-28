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
        real(rk),       intent(in)  :: node(3)

        !real(rk)    :: xi, eta, zeta
        real(rk)    :: polyval


        !xi   = node%c1_
        !eta  = node%c2_
        !zeta = node%c3_

        polyval = LegendreVal(spacedim,mode,node(1),node(2),node(3))

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
        real(rk),       intent(in)  :: node(3)
        integer(ik),    intent(in)  :: dir

        !real(rk)    :: xi, eta, zeta
        real(rk)    :: dpolyval

        !xi   = node%c1_
        !eta  = node%c2_
        !zeta = node%c3_

        dpolyval = DLegendreVal(spacedim,mode,node(1),node(2),node(3),dir)

    end function DPolynomialVal
    !****************************************************************************************





end module mod_polynomial
