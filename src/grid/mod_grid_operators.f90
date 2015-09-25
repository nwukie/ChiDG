module mod_grid_operators
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use type_point,         only: point_t
    use type_element,       only: element_t
    use type_domain,        only: domain_t
    use type_blockvector,   only: blockvector_t
    use type_densevector,   only: densevector_t
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

        integer(ik)             :: ielem, ierr
        real(rk), allocatable   :: fmodes(:)


        ! Check that variable index 'ivar' is valid
        if (ivar > domain%eqnset%neqns ) call signal(FATAL,'initialize_variable: variable index ivar exceeds the number of equations')


        !
        ! Loop through elements in mesh and call function projection
        !
        do ielem = 1,domain%mesh%nelem
            associate (elem  =>  domain%mesh%elems(ielem), q => domain%sdata%q%lvecs(ielem))

                ! Initial array allocation
                if (.not. allocated(fmodes)) allocate(fmodes(q%nterms()))


                ! Reallocate mode storage if necessary. For example, if the order of the expansion was changed
                if (size(fmodes) /= q%nterms()) then
                    if (allocated(fmodes)) deallocate(fmodes)
                    allocate(fmodes(q%nterms()), stat=ierr)
                    if (ierr /= 0) call AllocationError
                end if


                if (.not. allocated(fmodes)) call signal(FATAL,"initialize_variable: fmodes not allocated")



                !
                ! Call function projection
                !
                call project_function_xyz(fcn,elem%nterms_s,elem%coords,fmodes)



                !
                ! Store the projected modes to the solution expansion
                !
                !q%mat(:,ivar) = fmodes
                call q%setvar(ivar,fmodes)

            end associate
        end do ! ielem

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


        !
        ! Evaluate polynomial modes at node location
        !
        do iterm = 1,elem%nterms_c

            polyvals(iterm) = polynomialVal(3,elem%nterms_c,iterm,node)

        end do


        !
        ! Evaluate mesh point from dot product of modes and polynomial values
        !
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
        class(densevector_t), intent(inout)  :: q
        integer(ik),        intent(in)     :: ivar
        real(rk),           intent(in)     :: xi,eta,zeta

        real(rk)                   :: val
        type(point_t)              :: node
        real(rk)                   :: polyvals(q%nterms())
        integer(ik)                :: iterm, ielem


        call node%set(xi,eta,zeta)


        !
        ! Evaluate polynomial modes at node location
        !
        do iterm = 1,q%nterms()

            polyvals(iterm)  = polynomialVal(3,q%nterms(),iterm,node)

        end do


        !
        ! Evaluate x from dot product of modes and polynomial values
        !
        val = dot_product(q%getvar(ivar),polyvals)

    end function






end module mod_grid_operators
