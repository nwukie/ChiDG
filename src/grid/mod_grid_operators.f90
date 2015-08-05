module mod_grid_operators
    use mod_kinds,          only: rk, ik
    use type_point,         only: point_t
    use type_element,       only: element_t
    use type_domain,        only: domain_t
    use type_expansion,     only: expansion_t
    use atype_function,     only: function_t
    use mod_polynomial,     only: PolynomialVal
    use mod_project,        only: project_function_xyz
    implicit none


contains


    !>  Project functions to element solution variables
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  domain  Domain upon which the solution is projected
    !!  @param[in]  ivar    Integer index of the variable being initialized
    !!  @param[in]  fcn     Function being projected to the solution
    !---------------------------------------------------------------------------------
    subroutine initialize_variable(domain,ivar,fcn)
        type(domain_t),     intent(inout)   :: domain
        integer(ik),        intent(in)      :: ivar
        class(function_t),  intent(inout)   :: fcn

        integer(ik)             :: ielem
        real(rk), allocatable   :: fmodes(:)


        ! Loop through elements in mesh and call function projection
        do ielem = 1,domain%mesh%nelem
            associate (elem  =>  domain%mesh%elems(ielem), q => domain%sdata%q(ielem))

                ! Reallocate mode storage if necessary
                if (size(fmodes) /= q%nterms) then
                    if (allocated(fmodes)) deallocate(fmodes)
                    allocate(fmodes(q%nterms))
                end if

                ! Call function projection
                call project_function_xyz(fcn,elem%nterms_s,elem%coords,fmodes)

                ! Store the projected modes to the solution expansion
                q%mat(:,ivar) = fmodes

            end associate

        end do

    end subroutine







    !>  Compute a coordinate value, based on the location in reference space (xi, eta, zeta)
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  elem    Element containing coordinate expansion
    !!  @param[in]  icoord  Integer corresponding to coordinate index 1(x), 2(y), 3(z)
    !!  @param[in]  xi      Real value for xi-coordinate
    !!  @param[in]  eta     Real value for eta-coordinate
    !!  @param[in]  zeta    Real value for zeta-coordinate
    !----------------------------------------------------------------------------------
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





    !>  Compute a variable value, based on the location in reference space (xi, eta, zeta)
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  q       Solution expansion for a given element
    !!  @param[in]  ivar    Integer corresponding to variable index
    !!  @param[in]  xi      Real value for xi-coordinate
    !!  @param[in]  eta     Real value for eta-coordinate
    !!  @param[in]  zeta    Real value for zeta-coordinate
    !-----------------------------------------------------------------------------------
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
