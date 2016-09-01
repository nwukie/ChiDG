module mod_grid_operators
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: TWO_DIM, THREE_DIM, X_DIR, Y_DIR, Z_DIR, ZETA_DIR, ONE, ZERO
    use type_point,         only: point_t
    use type_element,       only: element_t
    use type_blockvector,   only: blockvector_t
    use type_densevector,   only: densevector_t
    use type_solverdata,    only: solverdata_t
    use type_chidg_data,    only: chidg_data_t
    use type_function,      only: function_t
    use mod_polynomial,     only: PolynomialVal, DPolynomialVal
    use mod_project,        only: project_function_xyz
    implicit none


contains


    !>  Project functions to element solution variables
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  domain  Domain upon which the solution is projected
    !!  @param[in]  ivar    Integer index of the variable being initialized
    !!  @param[in]  fcn     Function being projected to the solution
    !!
    !--------------------------------------------------------------------------------------------------------
    subroutine initialize_variable(data,ivar,fcn)
        type(chidg_data_t),     intent(inout)   :: data
        integer(ik),            intent(in)      :: ivar
        class(function_t),      intent(inout)   :: fcn

        integer(ik)             :: ielem, ierr, idom, nterms, spacedim
        real(rk), allocatable   :: fmodes(:)


        !
        ! Loop through elements in mesh and call function projection
        !
        do idom = 1,data%ndomains()
            ! Check that variable index 'ivar' is valid
            if (ivar > data%eqnset(idom)%prop%nequations() ) call chidg_signal(FATAL,'initialize_variable: variable index ivar exceeds the number of equations')

            do ielem = 1,data%mesh(idom)%nelem

                    !
                    ! Get spacedim
                    !
                    spacedim = data%mesh(idom)%elems(ielem)%spacedim

                    !
                    ! Initial array allocation
                    !
                    nterms = data%sdata%q%dom(idom)%vecs(ielem)%nterms()
                    if (.not. allocated(fmodes)) allocate(fmodes(nterms))


                    !
                    ! Reallocate mode storage if necessary. For example, if the order of the expansion was changed
                    !
                    if (size(fmodes) /= nterms) then
                        if (allocated(fmodes)) deallocate(fmodes)
                        allocate(fmodes(nterms), stat=ierr)
                        if (ierr /= 0) call AllocationError
                    end if


                    if (.not. allocated(fmodes)) call chidg_signal(FATAL,"initialize_variable: fmodes not allocated")


                    !
                    ! Call function projection
                    !
                    call project_function_xyz(fcn,spacedim,data%mesh(idom)%elems(ielem)%nterms_s,data%mesh(idom)%elems(ielem)%coords,fmodes)


                    !
                    ! Store the projected modes to the solution expansion
                    !
                    call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ivar,fmodes)

            end do ! ielem

        end do ! idomain

    end subroutine initialize_variable
    !**************************************************************************************************************












end module mod_grid_operators
