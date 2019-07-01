module type_harmonic_balance
#include<messenger.h>
    use messenger,                      only: write_line
    use mod_kinds,                      only: rk,ik
    use mod_constants,                  only: NO_FACE, NO_ID, ZERO
    use mod_spatial,                    only: update_space

    use type_time_integrator_spectral,  only: time_integrator_spectral_t
    use type_system_assembler,          only: system_assembler_t
    
    use type_chidg_data,                only: chidg_data_t
    use type_nonlinear_solver,          only: nonlinear_solver_t
    use type_linear_solver,             only: linear_solver_t
    use type_preconditioner,            only: preconditioner_t
    use type_rvector,                   only: rvector_t
    use mod_HB_matrices,                only: calc_pseudo_spectral_operator
    use type_seed,                      only: seed_t
    use type_element_info,              only: element_info_t, element_info
    use type_chidg_vector
    use DNAD_D

    implicit none
    private


    !>
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   2/13/2017
    !!
    !---------------------------------------------------------------------------------
    type, extends(time_integrator_spectral_t), public    :: harmonic_balance_t
    
    contains
        procedure   :: init
        procedure   :: step
    end type harmonic_balance_t
    !*********************************************************************************



    !>  Object for assembling the implicit system
    !!  
    !!  @author Matteo Ugolotti + Mayank Sharma
    !!  @date 2/13/2017
    !!
    !!
    !---------------------------------------------------------------------------------
    type , extends(system_assembler_t), public :: assemble_harmonic_balance_t

    contains
        procedure   :: assemble
    end type assemble_harmonic_balance_t
    !*********************************************************************************




contains



    !>  Initialize the harmonic_balance_t time integrator 
    !!  
    !!  Main activity here, crate the assembler and attach it to the time_integrator
    !!  object so it can be passed to the nonlinear solver.
    !!
    !!  @author Matteo Ugolotti + Mayank Sharma
    !!  @date 2/13/2017
    !!
    !!  @author Mayank Sharma
    !!  @date   2/26/2017
    !!
    !---------------------------------------------------------------------------------
    subroutine init(self)
        class(harmonic_balance_t),  intent(inout)   :: self

        integer(ik)                         :: ierr
        type(assemble_harmonic_balance_t)   :: assemble_harmonic_balance

        call self%set_name('Harmonic Balance')

        if (allocated(self%system)) deallocate(self%system)
        allocate(self%system, source=assemble_harmonic_balance, stat=ierr)
        if (ierr /= 0) call AllocationError

    end subroutine init
    !*********************************************************************************




    !>  Solving a steady-like equation set. Time levels are included in the LHS
    !!  by means of pseudo-spectral operator
    !!
    !!  Solving R(Q) = 0   
    !!
    !!  @author Matteo Ugolotti + Mayank Sharma
    !!  @date 2/13/2017
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine step(self,data,nonlinear_solver,linear_solver,preconditioner)
        class(harmonic_balance_t),                      intent(inout)   :: self
        type(chidg_data_t),                             intent(inout)   :: data
        class(nonlinear_solver_t),   optional,          intent(inout)   :: nonlinear_solver
        class(linear_solver_t),      optional,          intent(inout)   :: linear_solver
        class(preconditioner_t),     optional,          intent(inout)   :: preconditioner

        ! Simply solve the nonlinear system. No iteration in time.
        call nonlinear_solver%solve(data,self%system,linear_solver,preconditioner)

        ! Store end residual from nonlinear solver.
        call self%residual_norm%push_back(nonlinear_solver%residual_norm%at(nonlinear_solver%residual_norm%size()))

    end subroutine step
    !*********************************************************************************



    !>  Assemble the system for the 'Harmonic Balance' equations with temporal contributions
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti
    !!  @date   2/13/2017
    !!
    !---------------------------------------------------------------------------------
    subroutine assemble(self,data,differentiate,timing)
        class(assemble_harmonic_balance_t),   intent(inout)               :: self
        type(chidg_data_t),                   intent(inout)               :: data
        logical,                              intent(in)                  :: differentiate
        real(rk),                             intent(inout), optional     :: timing

        integer(ik)             :: itime_outer, itime_inner, idom, ielem, ivar, eqn_ID, ierr, ntime, nfields, &
                                   irow_start, irow_end, icol_start, icol_end, nterms_s
        real(rk),   allocatable :: temp_1(:), temp_2(:), D(:,:), temp_mat(:,:), hb_mat(:,:)
        type(chidg_vector_t)    :: rhs_tmp
        type(seed_t)            :: seed
        type(element_info_t)    :: elem_info

        
        associate ( rhs => data%sdata%rhs, lhs => data%sdata%lhs, q => data%sdata%q )
        call rhs%clear()
        if (differentiate) call lhs%clear()

        ! Set local variables equal to the values set in time_manager
        ntime = size(data%time_manager%times)
        D = data%time_manager%D

        do itime_outer = 1,ntime

            ! Spatial update needed
            data%time_manager%itime = itime_outer
            data%time_manager%t     = data%time_manager%times(itime_outer)
            call update_space(data,differentiate,timing)

            do idom = 1,data%mesh%ndomains()

                if (allocated(temp_1) .and. allocated(temp_2)) deallocate(temp_1,temp_2)
                allocate(temp_1(data%mesh%domain(idom)%nterms_s),temp_2(data%mesh%domain(idom)%nterms_s), stat=ierr)
                if (ierr /= 0) call AllocationError

                do ielem = 1,data%mesh%domain(idom)%nelem
                    eqn_ID = data%mesh%domain(idom)%elems(ielem)%eqn_ID
                    nterms_s = data%mesh%domain(idom)%elems(ielem)%nterms_s
                    nfields = data%eqnset(eqn_ID)%prop%nprimary_fields()
                    do itime_inner = 1,ntime

                        ! LHS HB contribution for each variable
                        temp_mat = D(itime_outer,itime_inner)*data%mesh%domain(idom)%elems(ielem)%mass
                        if (allocated(hb_mat)) deallocate(hb_mat)
                        allocate(hb_mat(nterms_s*nfields,nterms_s*nfields), stat=ierr)
                        if (ierr /= 0) call AllocationError
                        hb_mat(:,:) = ZERO

                        ! Accumulate variable contributions
                        do ivar = 1,data%eqnset(eqn_ID)%prop%nprimary_fields()

                            ! Temporary variables for computing the temporal rhs contributions
                            temp_1 = D(itime_outer,itime_inner)*matmul(data%mesh%domain(idom)%elems(ielem)%mass, &
                                                            q%dom(idom)%vecs(ielem)%getvar(ivar,itime_inner))

                            temp_2 = rhs%dom(idom)%vecs(ielem)%getvar(ivar,itime_outer) + temp_1

                            ! Store the temporal contributions in the rhs
                            call rhs%dom(idom)%vecs(ielem)%setvar(ivar,itime_outer,temp_2)

                            ! Accumulate hb contribution for the variable
                            irow_start = 1 + (ivar-1)*nterms_s
                            irow_end   = irow_start + (nterms_s-1)
                            icol_start = 1 + (ivar-1)*nterms_s
                            icol_end   = icol_start + (nterms_s-1)
                            hb_mat(irow_start:irow_end,icol_start:icol_end) = temp_mat

                        end do  ! ivar


                        ! LHS contribution
                        elem_info = data%mesh%get_element_info(idom,ielem)



                        call seed%init(data%mesh%domain(idom)%elems(ielem)%idomain_g,  &
                                       data%mesh%domain(idom)%elems(ielem)%idomain_l,  &
                                       data%mesh%domain(idom)%elems(ielem)%ielement_g, &
                                       data%mesh%domain(idom)%elems(ielem)%ielement_l, &
                                       data%mesh%domain(idom)%elems(ielem)%nfields,    &
                                       data%mesh%domain(idom)%elems(ielem)%nterms_s,   &
                                       IRANK,                                          &
                                       itime_inner,                                    &
                                       data%mesh%domain(idom)%elems(ielem)%dof_start,  &
                                       NO_ID,                                          &
                                       NO_ID,                                          &
                                       NO_ID)

                        ! Store HB contribution for all fields to lhs at one time
                        if (itime_inner /= itime_outer) then
                            call lhs%dom(idom)%store_hb_element(hb_mat,elem_info,seed,itime_outer)
                        end if


                    end do  ! itime_inner
                end do  ! ielem
            end do  ! idom
        end do  ! itime_outer

        end associate


    end subroutine assemble
    !*********************************************************************************






end module type_harmonic_balance
