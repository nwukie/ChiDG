module mod_grid_operators
    use mod_kinds,          only: rk, ik
    use type_point,         only: point_t
    use type_element,       only: element_t
    use type_expansion,     only: expansion_t
    use mod_polynomial,     only: PolynomialVal
    implicit none


contains




    function mesh_point(elem,icoord,xi,eta,zeta) result(val)
        class(element_t),   intent(in)  :: elem
        integer(ik),        intent(in)  :: icoord
        real(rk),           intent(in)  :: xi, eta, zeta

        real(rk)        :: val
        type(point_t)   :: node
        real(rk)        :: polyvals(elem%nterms_c)
        integer(ik)     :: iterm, ielem

        if (icoord > 3) stop "Error: mesh_point -- icoord exceeded 3 physical coordinates"

        call node%set(xi,eta,zeta)

        ! Evaluate polynomial modes at node location
        do iterm = 1,elem%nterms_c
            polyvals(iterm) = polynomialVal(3,elem%nterms_c,iterm,node)
        end do

        ! Evaluate mesh point from dot product of modes and polynomial values
        val = dot_product(elem%coords%mat(:,icoord), polyvals)

    end function




    function solution_point(q,ivar,xi,eta,zeta) result(val)
        class(expansion_t), intent(inout)  :: q
        integer(ik),        intent(in)     :: ivar
        real(rk),           intent(in)     :: xi,eta,zeta

        real(rk)                   :: val
        type(point_t)              :: node
        real(rk)                   :: polyvals(q%nterms)
        integer(ik)                :: iterm, ielem


        call node%set(xi,eta,zeta)

        ! Evaluate polynomial modes at node location
        do iterm = 1,q%nterms
            polyvals(iterm)  = polynomialVal(3,q%nterms,iterm,node)
        end do

        ! Evaluate x from dot product of modes and polynomial values
        val = dot_product(q%var(ivar),polyvals)
    end function






end module mod_grid_operators
