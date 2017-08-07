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
    function polynomial_val(spacedim,nterms,mode,node) result(polyval)
        integer(ik),    intent(in)  :: spacedim
        integer(ik),    intent(in)  :: nterms
        integer(ik),    intent(in)  :: mode
        real(rk),       intent(in)  :: node(3)

        real(rk)    :: polyval

        polyval = legendre_val(spacedim,mode,node(1),node(2),node(3))

    end function polynomial_val
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
    function dpolynomial_val(spacedim,nterms,mode,node,dir) result(dpolyval)
        integer(ik),    intent(in)  :: spacedim
        integer(ik),    intent(in)  :: nterms
        integer(ik),    intent(in)  :: mode
        real(rk),       intent(in)  :: node(3)
        integer(ik),    intent(in)  :: dir

        real(rk)    :: dpolyval

        dpolyval = dlegendre_val(spacedim,mode,node(1),node(2),node(3),dir)

    end function dpolynomial_val
    !****************************************************************************************





    !>  Return second/mixed derivative of a term in a polynomial expansion at a specified 
    !!  (xi, eta, zeta) location.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/7/2017
    !!
    !!  @param[in]  spacedim   Integer specifying whether to compute the 1D, 2D, 
    !!                          or 3D expansion function.
    !!  @param[in]  nterms      Number of terms in the polynomial expansion.
    !!  @param[in]  mode        Integer of the mode in the expansion to be computed.
    !!  @param[in]  node        point_t containing (xi,eta,zeta) location on reference element.
    !!  @param[in]  dir         Integer indicating direction of the derivative (xi, eta, zeta).
    !!
    !----------------------------------------------------------------------------------------
    function ddpolynomial_val(spacedim,nterms,mode,node,partial1,partial2) result(res)
        integer(ik),    intent(in)  :: spacedim
        integer(ik),    intent(in)  :: nterms
        integer(ik),    intent(in)  :: mode
        real(rk),       intent(in)  :: node(3)
        integer(ik),    intent(in)  :: partial1
        integer(ik),    intent(in)  :: partial2

        real(rk)    :: res

        res = ddlegendre_val(spacedim,mode,node(1),node(2),node(3),partial1,partial2)

    end function ddpolynomial_val
    !****************************************************************************************








end module mod_polynomial
