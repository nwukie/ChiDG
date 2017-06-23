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
    use type_chidg_vector
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
        integer(ik)                 :: niter
        real(rk)                    :: resid, resid_prev, timing, entropy_error, residual_ratio
        real(rk), allocatable       :: vals(:)
        type(chidg_vector_t)        :: b, qn, qold, qnew

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
            qn         = q      ! Store qn, since q will be operated on in the inner loop
            resid      = ONE    ! Force inner loop entry
            resid_prev = ONE    !
            niter      = 0      ! Initialize inner loop counter


            !
            ! NONLINEAR CONVERGENCE LOOP
            !
            do while ( resid > self%tol )
                niter = niter + 1


                !
                ! Store the value of the current inner iteration solution (k) for the solution update (n+1), q_(n+1)_k
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
                ! Tolerance check
                !
                call self%residual_norm%push_back(resid)
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
                ! Advance solution with update vector
                !
                q = qold + dq



                ! Clear working storage
                call dq%clear()


                ! Print iteration information
                call write_line(niter, resid, linear_solver%niter, controller%lhs_updated, controller%preconditioner_updated, delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))


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
