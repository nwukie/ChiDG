module type_newton
    use messenger,              only: write_line
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, DIAG
    use mod_chidg_mpi,          only: ChiDG_COMM, GLOBAL_MASTER
    use mod_io,                 only: verbosity
    use type_nonlinear_solver,  only: nonlinear_solver_t
    use type_linear_solver,     only: linear_solver_t
    use type_chidg_data,        only: chidg_data_t
    use type_system_assembler,  only: system_assembler_t
    use type_preconditioner,    only: preconditioner_t
    use type_solver_controller, only: solver_controller_t
    use mod_tecio_old,          only: write_tecio_old
    use type_chidg_vector
    use ieee_arithmetic,        only: ieee_is_nan
    implicit none
    private



    !>  Solution advancement via the newton's method
    !!
    !!
    !!  @author Nathan A. Wukie
    !!
    !----------------------------------------------------------------------------------------
    type, extends(nonlinear_solver_t), public :: newton_t


    contains

        procedure   :: solve

    end type newton_t
    !****************************************************************************************






contains





    !>  Solve for Newton update 'dq'
    !!
    !!  @author Nathan A. Wukie
    !!
    !-----------------------------------------------------------------------------------------
    subroutine solve(self,data,system,linear_solver,preconditioner,solver_controller)
        class(newton_t),                            intent(inout)           :: self
        type(chidg_data_t),                         intent(inout)           :: data
        class(system_assembler_t),      optional,   intent(inout)           :: system
        class(linear_solver_t),         optional,   intent(inout)           :: linear_solver
        class(preconditioner_t),        optional,   intent(inout)           :: preconditioner
        class(solver_controller_t),     optional,   intent(inout),  target  :: solver_controller

        character(100)              :: filename
        integer(ik)                 :: niter, step
        real(rk)                    :: resid, resid0, resid_prev, resid_new, alpha, f0, fn, forcing_term, &
                                       timing, entropy_error, residual_ratio
        logical                     :: absolute_convergence, relative_convergence, searching
        type(chidg_vector_t)        :: b, qn, qold, qnew, q0




        type(solver_controller_t),  target  :: default_controller
        class(solver_controller_t), pointer :: controller


        if (present(solver_controller)) then
            controller => solver_controller
        else
            controller => default_controller
        end if




        associate ( q   => data%sdata%q,    &
                    dq  => data%sdata%dq,   &
                    rhs => data%sdata%rhs,  &
                    lhs => data%sdata%lhs)

            call write_line('NONLINEAR SOLVER', io_proc=GLOBAL_MASTER, silence=(verbosity<2))
            call write_line("iter","|R(Q)|","Linear Solver(niter)", "LHS Updated", "Preconditioner Updated", delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))



            !
            ! start timer
            !
            call self%timer%start()

            
            !
            ! Startup values
            !
            absolute_convergence = .true.
            relative_convergence = .true.
            qn                   = q      ! Store qn, since q will be operated on in the inner loop
            resid                = ONE    ! Force inner loop entry
            resid0               = ONE
            resid_prev           = ONE    
            niter                = 0      ! Initialize inner loop counter


            !
            ! NONLINEAR CONVERGENCE LOOP
            !
            !do while ( resid > self%tol )
            do while( absolute_convergence .and. relative_convergence )
                niter = niter + 1


                !
                ! Store the value of the current inner iteration solution (k) for the solution update (n+1), q_(n+1)_k
                !
                qold = q


                !
                ! Update Spatial Residual and Linearization (rhs, lin)
                !
                if ( niter <= 2)  then
                    residual_ratio = ONE
                else
                    residual_ratio = resid/resid_prev
                end if

                call system%assemble(data, timing=timing, differentiate=controller%update_lhs(lhs,niter,residual_ratio) )


                !
                ! Compute residual norms
                !
                resid_prev = resid
                resid      = rhs%norm(ChiDG_COMM)
                if (niter == 1) resid0 = resid


                !
                ! Tolerance check
                !
                call self%residual_norm%push_back(resid)
                call write_line("|R| = ", resid, io_proc=GLOBAL_MASTER, silence=(verbosity<4))
                if ( resid < self%tol ) then
                    call write_line(niter, resid, 0, controller%lhs_updated, .false., delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))
                    exit
                end if
                call self%residual_time%push_back(timing)   ! non-essential record-keeping






                !
                ! Assign rhs to b, which should allocate storage
                !
                !b = (rhs)  ! BEWARE: this causes an error. Parentheses operator not defined
                b = (-ONE)*rhs




                !
                ! We need to solve the matrix system Ax=b for the update vector x (dq)
                !
                call linear_solver%solve(lhs,dq,b,preconditioner,controller)

                call self%matrix_iterations%push_back(linear_solver%niter)
                call self%matrix_time%push_back(linear_solver%timer%elapsed())




                !
                ! Line Search for appropriate step
                !
                q0 = qold
                f0 = resid
                step = 0
                if (self%search) then

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
                        call system%assemble(data,timing=timing,differentiate=.false.)

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
                q = qn
                call dq%clear()



                ! Print iteration information
                call write_line(niter, resid, linear_solver%niter, controller%lhs_updated, controller%preconditioner_updated, delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))



                !
                ! Convergence switches
                !
                absolute_convergence = (resid > self%tol)
                relative_convergence = ( (log10(resid0) - log10(resid)) < real(self%norders_reduction,rk) )


                call data%sdata%q_out%init(data%mesh,data%time_manager%ntime)
                call data%sdata%q_out%set_ntime(data%time_manager%ntime)
                call data%sdata%q_out%clear()
                data%sdata%q_out = data%sdata%q

                write(filename, "(A,F8.6,A4)") 'quasi_newton_', real(niter,rk), '.plt'
                call write_tecio_old(data,filename, write_domains=.true., write_surfaces=.false.)




            end do ! while error


            !
            ! stop timer. Record timings
            !
            call self%timer%stop()
            call write_line('Nonlinear solver elapsed time: ', self%timer%elapsed(), io_proc=GLOBAL_MASTER, silence=(verbosity<3))
            !call self%timer%report('Solver Elapsed Time:')
            call self%total_time%push_back(self%timer%elapsed())



        end associate



        !
        ! Store newton iteration count
        !
        call self%newton_iterations%push_back(niter)

    end subroutine solve
    !****************************************************************************************





end module type_newton
