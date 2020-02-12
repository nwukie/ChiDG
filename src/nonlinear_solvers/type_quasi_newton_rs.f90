module type_quasi_newton_rs
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, DIAG, NO_DIFF
    use mod_chidg_mpi,          only: ChiDG_COMM, GLOBAL_MASTER, IRANK, NRANK
    use mod_io,                 only: verbosity
    use mpi_f08,                only: MPI_Barrier

    use type_chidg_data,        only: chidg_data_t
    use type_nonlinear_solver,  only: nonlinear_solver_t
    use type_system_assembler,  only: system_assembler_t
    use type_linear_solver,     only: linear_solver_t
    use type_preconditioner,    only: preconditioner_t
    use type_solver_controller, only: solver_controller_t
    use type_chidg_vector

    use ieee_arithmetic,        only: ieee_is_nan
    implicit none
    private



    !>  Solution advancement via the backward-euler method
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !------------------------------------------------------------------------------------------
    type, extends(nonlinear_solver_t), public :: quasi_newton_rs_t


    contains
    
        procedure   :: solve
        final       :: destructor

    end type quasi_newton_rs_t
    !******************************************************************************************










contains





    !>  Solve for update 'dq'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine solve(self,data,system,linear_solver,preconditioner,solver_controller)
        class(quasi_newton_rs_t),               intent(inout)           :: self
        type(chidg_data_t),                     intent(inout)           :: data
        class(system_assembler_t),  optional,   intent(inout)           :: system
        class(linear_solver_t),     optional,   intent(inout)           :: linear_solver
        class(preconditioner_t),    optional,   intent(inout),  target  :: preconditioner
        class(solver_controller_t), optional,   intent(inout),  target  :: solver_controller

        type(chidg_data_t)                                              :: data_copy
        character(100)          :: filename
        integer(ik)             :: itime, nsteps, ielem, iblk, iindex,  &
                                   niter, ieqn, idom, ierr,                     &
                                   rstart, rend, cstart, cend, nterms, imat, iwrite, step, eqn_ID, icfl

        real(rk)                :: dtau, amp, cfl, timing, resid, resid_prev, resid0, resid_new,    &
                                   alpha, f0, fn, forcing_term, residual_ratio
        real(rk), allocatable   :: vals(:), cfln(:), rnorm0(:), rnorm(:)
        type(chidg_vector_t)    :: b, qn, qold, qnew, dqdtau, q0
        logical                 :: searching, absolute_convergence, relative_convergence

        type(solver_controller_t),  target  :: default_controller
        class(solver_controller_t), pointer :: controller
      

        if (present(solver_controller)) then
            controller => solver_controller
        else
            controller => default_controller
        end if


      
        call write_line('NONLINEAR SOLVER', io_proc=GLOBAL_MASTER, silence=(verbosity<2))
        call write_line("iter","|R(Q)|","CFL", "Linear Solver(niter)", "LHS Updated", "Preconditioner Updated", delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))


        associate ( q   => data%sdata%q,    &
                    dq  => data%sdata%dq,   &
                    rhs => data%sdata%rhs,  &
                    lhs => data%sdata%lhs)


            !
            ! start timer
            !
            call self%timer%reset()
            call self%timer%start()


            !
            ! Startup values
            !
            absolute_convergence = .true.
            relative_convergence = .true.
            qn     = q      ! Store qn, since q will be operated on in the inner loop
            resid  = ONE    ! Force inner loop entry
            niter  = 0      ! Initialize inner loop counter


            !
            ! NONLINEAR CONVERGENCE INNER LOOP
            !
            do while ( absolute_convergence .and. relative_convergence )
                niter = niter + 1


                !
                ! Store the value of the current inner iteration solution (k) 
                ! for the solution update (n+1), q_(n+1)_k
                !
                qold = q


                !
                ! Update Spatial Residual and Linearization (rhs, lin)
                !
                !if (present(solver_controller)) then
                !    if ( niter <= 2)  then
                !        residual_ratio = ONE
                !    else
                !        residual_ratio = resid/resid_prev
                !    end if

                !    call system%assemble( data,             &
                !                          timing=timing,    &
                !                          differentiate=solver_controller%update_lhs(niter,residual_ratio) )
                !else
                !    call system%assemble(data,timing=timing,differentiate=.true.)
                !end if
                if ( niter <= 2)  then
                    residual_ratio = ONE
                else
                    residual_ratio = resid/resid_prev
                end if

                call system%assemble( data,             &
                                      timing=timing,    &
                                      differentiate=controller%update_lhs(lhs,niter,residual_ratio) )
                resid_prev = resid
                resid      = rhs%norm(ChiDG_COMM)



                !
                ! Print diagnostics, check tolerance.
                !
                call self%residual_norm%push_back(resid)
                call write_line('|R| = ', resid, io_proc=GLOBAL_MASTER, silence=(verbosity<4))
                if ( resid < self%tol ) then
                    call write_line(niter, resid, cfln(1), 0, delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))
                    exit
                end if


                if ( ieee_is_nan(resid) ) then
                    call chidg_signal(FATAL,"quasi_newton_rs%solve: NaN residual calculation. Check initial solution and operator objects.")
                end if

                call self%residual_time%push_back(timing)


                !
                ! Compute and store residual norm for each field
                !
                if (niter == 1) then
                    resid0 = rhs%norm(ChiDG_COMM)
                    rnorm0 = rhs%norm_fields(ChiDG_COMM)
                    rnorm0 = resid0 !override
                end if

                rnorm = rhs%norm_fields(ChiDG_COMM)
                rnorm = resid !override




                !
                ! During the first few iterations, allow the initial residual norm to update
                ! if it has increased. Otherwise, if a solution converged to an error floor
                ! and the residual raised a little bit, then the CFL would essentially reset
                ! from infinity to CFL0, which we do not want.
                !
                if (niter < 5) then
                    where (rnorm > rnorm0)
                        rnorm0 = rnorm0 + 0.3_rk*(rnorm - rnorm0)
                    end where
                end if



                !
                ! Compute new cfl for each field
                !
                if (allocated(cfln)) deallocate(cfln)
                allocate(cfln(size(rnorm)), stat=ierr)
                if (ierr /= 0) call AllocationError

                where (rnorm /= 0.)
                    cfln = self%cfl0*(rnorm0/rnorm)
                else where
                    cfln = 0.1
                end where

                !if (IRANK == GLOBAL_MASTER) then
                !    call add_to_line("  CFL: ")
                !    do icfl = 1,size(cfln)
                !        call add_to_line(cfln(icfl))
                !    end do
                !    call send_line()
                !end if



                !
                ! Compute element-local pseudo-timestep
                !
                call compute_pseudo_timestep(data,cfln)


                !
                ! Add mass/dt to sub-block diagonal in dR/dQ
                !
                do idom = 1,data%mesh%ndomains()
                    do ielem = 1,data%mesh%domain(idom)%nelem
                        eqn_ID = data%mesh%domain(idom)%elems(ielem)%eqn_ID
                        do itime = 1,data%mesh%domain(idom)%ntime
                            do ieqn = 1,data%eqnset(eqn_ID)%prop%nprimary_fields()

                                ! get element-local timestep
                                dtau = data%mesh%domain(idom)%elems(ielem)%dtau(ieqn)

                                ! Need to compute row and column extends in diagonal so we can
                                ! selectively apply the mass matrix to the sub-block diagonal
                                nterms = data%mesh%domain(idom)%elems(ielem)%nterms_s
                                rstart = 1 + (ieqn-1) * nterms
                                rend   = (rstart-1) + nterms
                                cstart = rstart                 ! since it is square
                                cend   = rend                   ! since it is square

                                ! Add mass matrix divided by dt to the block diagonal
                                imat = lhs%dom(idom)%lblks(ielem,itime)%get_diagonal()
                                lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(rstart:rend,cstart:cend)  =  lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(rstart:rend,cstart:cend)  +  data%mesh%domain(idom)%elems(ielem)%mass*(ONE/dtau)

                            end do !ieqn
                        end do !itime
                    end do !ielem
                end do !idom


                data_copy%sdata%rhs = rhs
                data_copy%sdata%lhs = lhs

                !
                ! Assign rhs to b, which should allocate storage
                !
                !b = (rhs)  ! BEWARE: this causes an error. Parentheses operator not defined
                b = (-ONE)*rhs




                !
                ! We need to solve the matrix system Ax=b for the update vector x (dq)
                !

                ! Set forcing term. Converge 4 orders, or 1.e-8
                !forcing_term = resid/10000._rk
                !linear_solver%tol = max(1.e-8_rk, forcing_term)

                ! Solve system for newton step, dq
                call linear_solver%solve(data_copy%sdata%lhs,dq,b,preconditioner,controller)
                !call linear_solver%solve(lhs,dq,b,preconditioner,controller)
                call self%matrix_iterations%push_back(linear_solver%niter)
                call self%matrix_time%push_back(linear_solver%timer%elapsed())




                !
                ! Line Search for appropriate step
                !
                q0 = qold
                f0 = resid
                step = 0
                if (trim(self%search) == 'Backtrack') then

                    searching = .true.
                    do while (searching)

                        !
                        ! Set line increment via backtracking.
                        !   Try: 1, 0.5, 0.25... 2^-i
                        !
                        alpha = TWO**(-real(step,rk)) 
                        call write_line("       Testing newton direction with 'alpha' = ", alpha, io_proc=GLOBAL_MASTER, silence=(verbosity<3))


                        !
                        ! Advance solution along newton direction
                        !
                        qn = q0 + alpha*dq


                        !
                        ! Clear working vector
                        !
                        !call rhs%clear()


                        !
                        ! Set working solution. Test residual at (q). Do not differentiate
                        !
                        q = qn
                        call system%assemble(data,timing=timing,differentiate=NO_DIFF)

                        !
                        ! Compute new function value
                        !
                        fn = rhs%norm(ChiDG_COMM)


                        !
                        ! Test for |R| increasing too much or NaN. 
                        !   If residual is reasonably large, still allow some growth.
                        !   If the residual is small enough, we don't want any growth.
                        ! 
                        !
                        if (ieee_is_nan(fn)) then
                            searching = .true.
                        else if ( (fn > 1.e-3_rk) .and. (fn > 2.0_rk*f0) ) then
                            searching = .true.
                        else if ( (fn < 1.e-3_rk) .and. (fn > f0) ) then
                            searching = .true.
                        else
                            searching = .false.
                        end if

                        call write_line("       Rn(Q) = ", fn, io_proc=GLOBAL_MASTER, silence=(verbosity<3))

                        step = step + 1

                    end do

                else
                    qn = q0 + dq
                end if


                !
                ! Accept new solution, clear working storage, iterate
                !
                dq = qn - q0


                rhs = data_copy%sdata%rhs 
                lhs = data_copy%sdata%lhs 
!
                ! Add mass/dt to sub-block diagonal in dR/dQ
                !
                do idom = 1,data%mesh%ndomains()
                    do ielem = 1,data%mesh%domain(idom)%nelem
                        eqn_ID = data%mesh%domain(idom)%elems(ielem)%eqn_ID
                        do itime = 1,data%mesh%domain(idom)%ntime
                            do ieqn = 1,data%eqnset(eqn_ID)%prop%nprimary_fields()

                                ! get element-local timestep
                                dtau = data%mesh%domain(idom)%elems(ielem)%dtau(ieqn)

                                ! Need to compute row and column extends in diagonal so we can
                                ! selectively apply the mass matrix to the sub-block diagonal
                                nterms = data%mesh%domain(idom)%elems(ielem)%nterms_s
                                rstart = 1 + (ieqn-1) * nterms
                                rend   = (rstart-1) + nterms
                                cstart = rstart                 ! since it is square
                                cend   = rend                   ! since it is square


                                ! Add mass matrix divided by dt to the block diagonal

                                imat = lhs%dom(idom)%lblks(ielem,itime)%get_diagonal()
                                rhs%dom(idom)%vecs(ielem)%vec(rstart:rend)  = rhs%dom(idom)%vecs(ielem)%vec(rstart:rend) -  matmul(data%mesh%domain(idom)%elems(ielem)%mass*(ONE/dtau),dq%dom(idom)%vecs(ielem)%vec(rstart:rend))
                                !call rhs%dom(idom)%vecs(ielem)%setvar(ieqn, itime,  &
                                !rhs%dom(idom)%vecs(ielem)%getvar(ieqn, itime) &
                                !- matmul(data%mesh%domain(idom)%elems(ielem)%mass*(ONE/dtau),dq%dom(idom)%vecs(ielem)%getvar(ieqn, itime)))

                            end do !ieqn
                        end do !itime
                    end do !ielem
                end do !idom

!
                ! Assign rhs to b, which should allocate storage
                !
                !b = (rhs)  ! BEWARE: this causes an error. Parentheses operator not defined
                b = (-ONE)*rhs




                !
                ! We need to solve the matrix system Ax=b for the update vector x (dq)
                !

                ! Set forcing term. Converge 4 orders, or 1.e-8
                !forcing_term = resid/10000._rk
                !linear_solver%tol = max(1.e-8_rk, forcing_term)

                ! Solve system for newton step, dq
                call linear_solver%solve(lhs,dq,b,preconditioner,controller)
                call self%matrix_iterations%push_back(linear_solver%niter)
                call self%matrix_time%push_back(linear_solver%timer%elapsed())

                !
                ! Line Search for appropriate step
                !
                q0 = qold
                f0 = resid
                step = 0
                if (trim(self%search) == 'Backtrack') then

                    searching = .true.
                    do while (searching)

                        !
                        ! Set line increment via backtracking.
                        !   Try: 1, 0.5, 0.25... 2^-i
                        !
                        alpha = TWO**(-real(step,rk)) 
                        call write_line("       Testing newton direction with 'alpha' = ", alpha, io_proc=GLOBAL_MASTER, silence=(verbosity<3))


                        !
                        ! Advance solution along newton direction
                        !
                        qn = q0 + alpha*dq


                        !
                        ! Clear working vector
                        !
                        !call rhs%clear()


                        !
                        ! Set working solution. Test residual at (q). Do not differentiate
                        !
                        q = qn
                        call system%assemble(data,timing=timing,differentiate=NO_DIFF)

                        !
                        ! Compute new function value
                        !
                        fn = rhs%norm(ChiDG_COMM)


                        !
                        ! Test for |R| increasing too much or NaN. 
                        !   If residual is reasonably large, still allow some growth.
                        !   If the residual is small enough, we don't want any growth.
                        ! 
                        !
                        if (ieee_is_nan(fn)) then
                            searching = .true.
                        else if ( (fn > 1.e-3_rk) .and. (fn > 2.0_rk*f0) ) then
                            searching = .true.
                        else if ( (fn < 1.e-3_rk) .and. (fn > f0) ) then
                            searching = .true.
                        else
                            searching = .false.
                        end if

                        call write_line("       Rn(Q) = ", fn, io_proc=GLOBAL_MASTER, silence=(verbosity<3))

                        step = step + 1

                    end do

                else
                    qn = q0 + dq
                end if





                q = qn
                call dq%clear()







                absolute_convergence = (resid > self%tol)
                relative_convergence = (fn/resid0 > self%rtol)
                !relative_convergence = ( (log10(resid0) - log10(resid)) < real(self%norders_reduction,rk) )



                ! Print iteration information
                call write_line(niter, resid, cfln(1), linear_solver%niter, controller%lhs_updated, controller%preconditioner_updated, delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))


                call MPI_Barrier(ChiDG_COMM,ierr)
            end do ! niter


            !
            ! stop timer
            !
            call self%timer%stop()
            call self%total_time%push_back(self%timer%elapsed())
            call write_line('Nonlinear Solver elapsed time: ', self%timer%elapsed(), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<3))



        end associate



        call self%newton_iterations%push_back(niter)


    end subroutine solve
    !*****************************************************************************************











    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/17/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute_pseudo_timestep(data,cfln)
        type(chidg_data_t),     intent(inout)   :: data
        real(rk),               intent(in)      :: cfln(:)

        integer(ik) :: idom, eqn_ID

        !
        ! Loop through elements and compute time-step function
        !
        do idom = 1,data%mesh%ndomains()

            eqn_ID = data%mesh%domain(idom)%elems(1)%eqn_ID
            call data%eqnset(eqn_ID)%compute_pseudo_timestep(idom,data%mesh,data%sdata,cfln,itime = 1)

        end do !idom

    end subroutine compute_pseudo_timestep
    !*****************************************************************************************










    
    subroutine destructor(self)
        type(quasi_newton_rs_t),      intent(in) :: self

    end subroutine




end module type_quasi_newton_rs






