module type_newton
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO
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



    !>  Nonlinear iteration using newton's method.
    !!
    !!  Optional algorithms:
    !!      : Pseudo-transient continuation (&nonlinear_solve   ptc=.true. /)
    !!      : Residual Smoothing            (&nonlinear_solve   smooth=.true. /)
    !!      : Backtracking                  (&nonlinear_solve   search='Backtrack' /)
    !!
    !!  Termination:
    !!      After each step, the algorithm searches the working directory for a file 'STOP'.
    !!      Using this approach, a user can terminate iteration by creating a file called STOP
    !!      in the working directory of the problem. (i.e. touch STOP)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !------------------------------------------------------------------------------------------
    type, extends(nonlinear_solver_t), public :: newton_t

    contains
        procedure   :: solve
        procedure   :: record_and_report
        procedure   :: update_cfl
    end type newton_t
    !******************************************************************************************


contains


    !>  Solve for update 'dq'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine solve(self,data,system,linear_solver,preconditioner,solver_controller)
        class(newton_t),                        intent(inout)           :: self
        type(chidg_data_t),                     intent(inout)           :: data
        class(system_assembler_t),  optional,   intent(inout)           :: system
        class(linear_solver_t),     optional,   intent(inout)           :: linear_solver
        class(preconditioner_t),    optional,   intent(inout)           :: preconditioner
        class(solver_controller_t), optional,   intent(inout),  target  :: solver_controller

        character(:),   allocatable :: user_msg
        integer(ik)             :: itime, niter, ierr, icfl
        real(rk)                :: cfl, timing, resid, resid_prev, resid0, resid_new,    &
                                   alpha, f0, fn, forcing_term, residual_ratio
        real(rk), allocatable   :: cfln(:), rnorm0(:), rnorm(:)
        type(chidg_vector_t)    :: b, qn, qold, q0, f_smooth
        logical                 :: absolute_convergence, relative_convergence, stop_run, iteration_convergence

        type(solver_controller_t),  target  :: default_controller
        class(solver_controller_t), pointer :: controller
      

        ! Default controller
        controller => default_controller
        if (present(solver_controller)) controller => solver_controller

      
        call write_line('NONLINEAR SOLVER', io_proc=GLOBAL_MASTER, silence=(verbosity<2))
        call write_line("iter","|R(Q)|","CFL", "Linear Solver(niter)", "LHS Updated", "Preconditioner Updated", &
                        delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))


        ! start timer
        call self%timer%reset()
        call self%timer%start()


        associate ( q   => data%sdata%q,    &
                    dq  => data%sdata%dq,   &
                    rhs => data%sdata%rhs,  &
                    lhs => data%sdata%lhs)


        ! Startup values
        absolute_convergence = .true.   ! measures absolute convergence
        relative_convergence = .true.   ! measures relative convergence
        iteration_convergence= .true.   ! enforces maximum number of iterations
        stop_run             = .false.  ! watches for file in directory so user can stop nicely with output
        qn     = q                      ! Store qn, since q will be operated on in the inner loop
        resid  = ONE                    ! Force inner loop entry
        niter  = 0                      ! Initialize inner loop counter


        ! Nonlinear iteration, while not converged
        do while ( absolute_convergence .and. &
                   relative_convergence .and. &
                   iteration_convergence .and. &
                   (.not. stop_run) )
            niter = niter + 1


            ! Store the value of the current inner iteration solution (k) 
            ! for the solution update (n+1), q_(n+1)_k
            qold = q


            ! Update Spatial Residual and Linearization (rhs, lin)
            residual_ratio = resid/resid_prev
            if ( niter <= 2) residual_ratio = ONE
            call system%assemble( data,             &
                                  timing=timing,    &
                                  differentiate=controller%update_lhs(lhs,niter,residual_ratio) )

            ! Assign rhs to b, which should allocate storage
            !b = (rhs)  ! BEWARE: this causes an error. Parentheses operator not defined
            b = (-ONE)*rhs


            ! Compute and store residual norm for each field
            resid_prev = resid
            resid      = rhs%norm(ChiDG_COMM)
            if (niter == 1) then
                resid0 = rhs%norm(ChiDG_COMM)
                rnorm0 = rhs%norm_fields(ChiDG_COMM)
                rnorm0 = resid0 !override
            end if
            rnorm = rhs%norm_fields(ChiDG_COMM)
            rnorm = resid !override


            ! During the first few iterations, allow the initial residual norm to update
            ! if it has increased. Otherwise, if a solution converged to an error floor
            ! and the residual raised a little bit, then the CFL would essentially reset
            ! from infinity to CFL0, which we do not want.
            if (niter < 5) then
                where (rnorm > rnorm0)
                    rnorm0 = rnorm0 + 0.3_rk*(rnorm - rnorm0)
                end where
            end if


            ! Compute new cfl for each field
            call self%update_cfl(rnorm0,rnorm,cfln)


            ! Convergence check
            call self%record_and_report(resid,timing,niter,cfln(1))
            if (resid < self%tol) exit

            if ( ieee_is_nan(resid) ) then
                user_msg = "newton%solve: NaN residual norm. Check initial solution and operator objects."
                call chidg_signal(FATAL,user_msg)
            end if


            ! Pseudo-transient continuation
            if (self%ptc) then
                call contribute_pseudo_temporal(data,cfln)
            end if ! self%ptc

            ! Residual-smoothing
            if (self%smooth) then
                f_smooth = compute_smoothed_residual(data,preconditioner,cfln,controller)
                b = b - f_smooth
            end if ! self%smooth


            ! Solve system [lhs][dq] = [b] for newton step: [dq]
            call set_forcing_terms(linear_solver)
            call linear_solver%solve(lhs,dq,b,preconditioner,controller,data)


            ! Line Search for appropriate step
            q0 = qold
            f0 = resid
            select case (trim(self%search))
                case('Backtrack')
                    call backtracking(data,system,q0,qn,f0,fn)
                case('none','')
                    ! Update state, update residual, compute residual norm 
                    qn = q0 + dq
                    data%sdata%q = qn
                    call system%assemble(data,differentiate=.false.)

                    rhs = rhs + f_smooth
                    fn = rhs%norm(ChiDG_COMM)

                case default
                    user_msg = "Invalid nonlinear search algorithm in newton iteration routine. Valid inputs: 'Backtrack', 'none'."
                    call chidg_signal_one(OOPS, user_msg, trim(self%search))
            end select


            ! Accept new solution, update cfl using new residual, clear working storage
            q = qn
            call self%update_cfl(rnorm0,[fn],cfln)
            call dq%clear()


            ! Record iteration data
            call self%matrix_iterations%push_back(linear_solver%niter)
            call self%matrix_time%push_back(linear_solver%timer%elapsed())
            call self%residual_norm%push_back(fn)
            call self%residual_time%push_back(timing)
            absolute_convergence  = (fn > self%tol)
            relative_convergence  = (fn/resid0 > self%rtol)
            iteration_convergence = (niter < self%nmax) .or. self%nmax <= 0
            inquire(file='STOP', exist=stop_run)


            ! Print iteration information
            call write_line(niter+1, fn, cfln(1), linear_solver%niter, controller%lhs_updated, controller%preconditioner_updated, delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))
            call MPI_Barrier(ChiDG_COMM,ierr)


        end do ! niter


        ! stop timer
        call self%timer%stop()
        call self%total_time%push_back(self%timer%elapsed())
        call self%newton_iterations%push_back(niter)
        call write_line('Nonlinear Solver elapsed time: ', self%timer%elapsed(), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<3))

        end associate


    end subroutine solve
    !*****************************************************************************************







    !-----------------------------------------------------------------------------------------
    !
    !   Supplemental and helper subroutines
    !       : compute_pseudo_timestep
    !       : contribute_pseudo_temporal
    !       : backtracking
    !       : set_forcing_terms
    !
    !-----------------------------------------------------------------------------------------



    !>  Update cfl for each equation that is used to compute the pseudo timestep for
    !!  pseudo-transient continuation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/17/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine update_cfl(self,rnorm0,rnorm,cfln)
        class(newton_t),            intent(in)      :: self
        real(rk),                   intent(in)      :: rnorm0(:)
        real(rk),                   intent(in)      :: rnorm(:)
        real(rk),   allocatable,    intent(inout)   :: cfln(:)

        integer(ik) :: icfl, ierr

        ! Compute new cfl for each field
        if (allocated(cfln)) deallocate(cfln)
        allocate(cfln(size(rnorm0)), stat=ierr)
        if (ierr /= 0) call AllocationError

        if (size(rnorm)==1) then
            ! rnorm might be coming from a backtracking algorithm and 
            ! only give us a scalar value for the entire vector.
            where (rnorm /= 0.)
                cfln = self%cfl0*(rnorm0/rnorm(1))
            else where
                cfln = 0.1
            end where
        else

            where (rnorm /= 0.)
                cfln = self%cfl0*(rnorm0/rnorm)
            else where
                cfln = 0.1
            end where
        end if

        do icfl = 1,size(cfln)
            if (self%cflmax > 0. .and. cfln(icfl) > self%cflmax) cfln(icfl) = self%cflmax
        end do


    end subroutine update_cfl
    !*****************************************************************************************



    !>  Compute pseudo timestep for pseudo-transient continuation algorithm
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/17/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute_pseudo_timestep(data,cfln)
        type(chidg_data_t),     intent(inout)   :: data
        real(rk),               intent(in)      :: cfln(:)

        integer(ik) :: idom, eqn_ID

        ! Loop through elements and compute time-step function
        do idom = 1,data%mesh%ndomains()
            eqn_ID = data%mesh%domain(idom)%elems(1)%eqn_ID
            call data%eqnset(eqn_ID)%compute_pseudo_timestep(idom,data%mesh,data%sdata,cfln,itime = 1)
        end do !idom

    end subroutine compute_pseudo_timestep
    !*****************************************************************************************




    !>  Contribute pseudo-transient continuation terms to linear system matrix.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/17/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine contribute_pseudo_temporal(data,cfln)
        type(chidg_data_t), intent(inout)   :: data
        real(rk),           intent(in)      :: cfln(:)

        integer(ik) :: idom, ielem, eqn_ID, itime, ieqn, nterms, rstart, rend, cstart, cend, imat
        real(rk)    :: dtau

        ! Compute element-local pseudo-timestep
        call compute_pseudo_timestep(data,cfln)

        ! Add mass/dt to sub-block diagonal in dR/dQ
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
                        imat = data%sdata%lhs%dom(idom)%lblks(ielem,itime)%get_diagonal()
                        data%sdata%lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(rstart:rend,cstart:cend)  =  data%sdata%lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(rstart:rend,cstart:cend)  +  data%mesh%domain(idom)%elems(ielem)%mass*(ONE/dtau)

                    end do !ieqn
                end do !itime
            end do !ielem
        end do !idom

    end subroutine contribute_pseudo_temporal
    !*************************************************************************************




    !>  Compute smoothed residual
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   01/15/2019
    !!
    !-----------------------------------------------------------------------------------------
    function compute_smoothed_residual(data,M,cfln,controller) result(f_smooth)
        type(chidg_data_t),         intent(inout)   :: data
        class(preconditioner_t),    intent(inout)   :: M
        real(rk),                   intent(in)      :: cfln(:)
        class(solver_controller_t), intent(inout)   :: controller

        integer(ik)             :: idom, ielem, eqn_ID, itime, ifield, i
        real(rk)                :: dtau
        real(rk),   allocatable :: field(:)
        type(chidg_vector_t)    :: f_smooth

        ! Compute element-local pseudo-timestep
        call compute_pseudo_timestep(data,cfln)


        ! Update smoother(preconditioner)
        if (controller%update_preconditioner(data%sdata%lhs)) call M%update(data%sdata%lhs,data%sdata%rhs)


        ! Apply smoothing
        f_smooth = data%sdata%rhs
        f_smooth = M%apply(data%sdata%lhs,f_smooth)
        
        ! I don't think this iteration procedure is correct
        ! for more than one iteration
        !do i = 1,1
        !    f_smooth = M%apply(data%sdata%lhs,f_smooth)
        !end do ! smooth


        ! Scale vector by (M/dtau)
        do idom = 1,data%mesh%ndomains()
            do ielem = 1,data%mesh%domain(idom)%nelem
                eqn_ID = data%mesh%domain(idom)%elems(ielem)%eqn_ID
                do itime = 1,data%mesh%domain(idom)%ntime
                    do ifield = 1,data%eqnset(eqn_ID)%prop%nprimary_fields()

                        ! get element-local timestep
                        dtau = data%mesh%domain(idom)%elems(ielem)%dtau(ifield)

                        ! Retrieve field
                        field = f_smooth%dom(idom)%vecs(ielem)%getvar(ifield,itime)

                        ! Scale field by (M/dtau)
                        field = matmul(data%mesh%domain(idom)%elems(ielem)%mass/dtau,field)

                        ! Store scaled field 
                        call f_smooth%dom(idom)%vecs(ielem)%setvar(ifield,itime,field)

                    end do !ifield
                end do !itime
            end do !ielem
        end do !idom

    end function compute_smoothed_residual
    !*************************************************************************************


























    !>  Backtracking search procedure.
    !!
    !!  If the newton update(dq) causes the residual to increase, scale the update
    !!  to be smaller until the residual achieves some reasonable value.
    !!
    !!  Scaling based on recursive bisection:
    !!     Try: 1, 0.5, 0.25... 2^-i
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/20/2017
    !!
    !-------------------------------------------------------------------------------------
    subroutine backtracking(data,system,q0,qn,f0,fn)
        type(chidg_data_t),         intent(inout)   :: data
        class(system_assembler_t),  intent(inout)   :: system
        type(chidg_vector_t),       intent(inout)   :: q0
        type(chidg_vector_t),       intent(inout)   :: qn
        real(rk),                   intent(inout)   :: f0
        real(rk),                   intent(inout)   :: fn

        real(rk)                :: alpha
        logical                 :: searching
        integer(ik)             :: step

        step = 0
        searching = .true.
        do while (searching)

            alpha = TWO**(-real(step,rk)) 
            call write_line("       Testing newton direction with 'alpha' = ", alpha, io_proc=GLOBAL_MASTER, silence=(verbosity<3))

            ! Advance solution along newton direction
            qn = q0 + alpha*data%sdata%dq

            ! Set working solution. Test residual at (q). Do not differentiate
            data%sdata%q = qn

            ! Compute new function value and norm
            call system%assemble(data,differentiate=.false.)
            fn = data%sdata%rhs%norm(ChiDG_COMM)

            ! Test for |R| increasing too much or NaN. 
            !   If residual is reasonably large, still allow some growth.
            !   If the residual is small enough, we don't want any growth.
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


    end subroutine backtracking
    !************************************************************************************






    !>
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine set_forcing_terms(linear_solver)
        class(linear_solver_t), intent(inout)   :: linear_solver

        ! Set forcing term. Converge 4 orders, or 1.e-8
        !forcing_term = resid/10000._rk
        !linear_solver%tol = max(1.e-8_rk, forcing_term)

    end subroutine set_forcing_terms
    !***********************************************************************************




    !>
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine record_and_report(self,resid,timing,niter,cfln)
        class(newton_t),    intent(inout)   :: self
        real(rk),           intent(in)      :: resid
        real(rk),           intent(in)      :: timing
        integer(ik),        intent(in)      :: niter
        real(rk),           intent(in)      :: cfln

        call write_line('|R| = ', resid, io_proc=GLOBAL_MASTER, silence=(verbosity<4))
        if ( resid < self%tol ) then
!            call self%timer%stop()
!            call self%total_time%push_back(self%timer%elapsed())
            call self%residual_norm%push_back(resid)
            call self%residual_time%push_back(timing)
            call self%matrix_iterations%push_back(0)
            call self%matrix_time%push_back(0._rk)
            call self%newton_iterations%push_back(niter)
            call write_line(niter, resid, cfln, 0, delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))
!            exit
        !end if
        ! Compute and store first residual norm 
        else if (niter == 1) then
            call self%residual_norm%push_back(resid)
            call self%residual_time%push_back(timing)
            call self%matrix_iterations%push_back(0)
            call self%matrix_time%push_back(0._rk)
            call write_line(niter, resid, cfln, 0, delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))
        end if

    end subroutine record_and_report
    !***********************************************************************************






end module type_newton






