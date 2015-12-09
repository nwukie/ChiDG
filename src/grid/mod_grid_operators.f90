module mod_grid_operators
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use type_point,         only: point_t
    use type_element,       only: element_t
    use type_blockvector,   only: blockvector_t
    use type_densevector,   only: densevector_t
    use type_solverdata,    only: solverdata_t
    use type_chidg_data,    only: chidg_data_t
    use atype_function,     only: function_t
    use mod_polynomial,     only: PolynomialVal, DPolynomialVal
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
    subroutine initialize_variable(data,ivar,fcn)
        type(chidg_data_t),     intent(inout)   :: data
        integer(ik),            intent(in)      :: ivar
        class(function_t),      intent(inout)   :: fcn

        integer(ik)             :: ielem, ierr, idom, nterms
        real(rk), allocatable   :: fmodes(:)


        !
        ! Loop through elements in mesh and call function projection
        !
        do idom = 1,data%ndomains()
            ! Check that variable index 'ivar' is valid
            if (ivar > data%eqnset(idom)%item%neqns ) call chidg_signal(FATAL,'initialize_variable: variable index ivar exceeds the number of equations')

            do ielem = 1,data%mesh(idom)%nelem
                !associate (elem => data%mesh(idom)%elems(ielem), q => data%sdata%q%dom(idom)%lvecs(ielem))

                    ! Initial array allocation
                    nterms = data%sdata%q%dom(idom)%lvecs(ielem)%nterms()
                    if (.not. allocated(fmodes)) allocate(fmodes(nterms))


                    ! Reallocate mode storage if necessary. For example, if the order of the expansion was changed
                    if (size(fmodes) /= nterms) then
                        if (allocated(fmodes)) deallocate(fmodes)
                        allocate(fmodes(nterms), stat=ierr)
                        if (ierr /= 0) call AllocationError
                    end if


                    if (.not. allocated(fmodes)) call chidg_signal(FATAL,"initialize_variable: fmodes not allocated")



                    !
                    ! Call function projection
                    !
                    !call project_function_xyz(fcn,elem%nterms_s,elem%coords,fmodes)
                    call project_function_xyz(fcn,data%mesh(idom)%elems(ielem)%nterms_s,data%mesh(idom)%elems(ielem)%coords,fmodes)



                    !
                    ! Store the projected modes to the solution expansion
                    !
                    !call q%setvar(ivar,fmodes)
                    call data%sdata%q%dom(idom)%lvecs(ielem)%setvar(ivar,fmodes)

                !end associate
            end do ! ielem

        end do ! idomain

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

        if (icoord > 3) call chidg_signal(FATAL,"Error: mesh_point -- icoord exceeded 3 physical coordinates")

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
        val = dot_product(elem%coords%getvar(icoord), polyvals)

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
        class(densevector_t),   intent(in)     :: q
        integer(ik),            intent(in)     :: ivar
        real(rk),               intent(in)     :: xi,eta,zeta

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









    !> Compute coordinate metric term at a given point in computational space
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  elem        element_t containing the geometry definition and data
    !!  @param[in]  cart_dir    Cartesian coordinate being differentiated
    !!  @param[in]  comp_dir    Computational coordinate being differentiated with respect to
    !!  @param[in]  xi          Computational coordinate - xi
    !!  @param[in]  eta         Computational coordinate - eta
    !!  @param[in]  zeta        Computational coordinate - zeta
    !!
    !-------------------------------------------------------------------
    function metric_point(elem,cart_dir,comp_dir,xi,eta,zeta) result(val)
        class(element_t),   intent(in)  :: elem
        integer(ik),        intent(in)  :: cart_dir
        integer(ik),        intent(in)  :: comp_dir
        real(rk),           intent(in)  :: xi, eta, zeta
        


        real(rk)        :: val
        type(point_t)   :: node
        real(rk)        :: polyvals(elem%nterms_c)
        integer(ik)     :: iterm, ielem

        if (cart_dir > 3) call chidg_signal(FATAL,"Error: mesh_point -- card_dir exceeded 3 physical coordinates")
        if (comp_dir > 3) call chidg_signal(FATAL,"Error: mesh_point -- comp_dir exceeded 3 physical coordinates")

        call node%set(xi,eta,zeta)


        !
        ! Evaluate polynomial modes at node location
        !
        do iterm = 1,elem%nterms_c

            polyvals(iterm) = dpolynomialVal(3,elem%nterms_c,iterm,node,comp_dir)

        end do


        !
        ! Evaluate mesh point from dot product of modes and polynomial values
        !
        val = dot_product(elem%coords%getvar(cart_dir), polyvals)



    end function












end module mod_grid_operators
