module type_newton
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, NO_DIFF
    use mod_chidg_mpi,          only: ChiDG_COMM, GLOBAL_MASTER, IRANK, NRANK
    use mod_io,                 only: verbosity
    use mpi_f08,                only: MPI_Barrier

    use type_chidg_data,        only: chidg_data_t
    use type_nonlinear_solver,  only: nonlinear_solver_t
    use type_system_assembler,  only: system_assembler_t
    use type_linear_solver,     only: linear_solver_t
    use type_preconditioner,    only: preconditioner_t
    use type_solver_controller, only: solver_controller_t
    use type_timer,             only: timer_t
    use type_element_info,      only: element_info_t, element_info
    use type_chidg_vector
    use operator_chidg_mv

    use precon_jacobi,          only: precon_jacobi_t

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
        procedure   :: backtracking
        procedure   :: record_and_report
        procedure   :: apply_residual_smoother
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
        class(preconditioner_t),    optional,   intent(inout),  target  :: preconditioner
        class(solver_controller_t), optional,   intent(inout),  target  :: solver_controller

        character(:),   allocatable :: user_msg
        integer(ik)             :: niter, ierr
        real(rk)                :: cfl, timing, resid, resid_prev, resid0, resid_new,    &
                                   alpha, f0, fn, forcing_term, residual_ratio
        real(rk), allocatable   :: cfln(:), fn_fields(:)
        type(chidg_vector_t)    :: b, qn, qold, q0, f_smooth
        logical                 :: absolute_convergence, relative_convergence, stop_run, iteration_convergence

        type(solver_controller_t),  target  :: default_controller
        class(solver_controller_t), pointer :: controller
        class(preconditioner_t),    pointer :: smoother => null()
        type(precon_jacobi_t),      target  :: jacobi
        type(timer_t)                       :: timer_linear


        ! Default controller
        controller => default_controller
        if (present(solver_controller)) controller => solver_controller

      
        call write_line('NONLINEAR SOLVER', io_proc=GLOBAL_MASTER, silence=(verbosity<2))
        call write_line("iter","|R(Q)|","CFL", "Linear Solver(niter)", "LHS Updated", "Preconditioner Updated", &
                        delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))


        ! start timer
        call self%timer%reset()
        call self%timer%start()
        call timer_linear%reset()


        ! Initialize smoother
        select case (self%smoother)
            case('preconditioner','Preconditioner')
                smoother => preconditioner
            case('jacobi','Jacobi')
                smoother => jacobi
                call smoother%init(data)
            case default
                call chidg_signal(FATAL,"newton%solve: invalid smoother. 'default' or 'jacobi'.")
        end select


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

            if (niter == 1) then
                resid0 = rhs%norm(ChiDG_COMM)
                cfln = rhs%norm_fields(ChiDG_COMM) !initialize
                cfln = self%cfl0
            end if
            resid_prev = resid
            resid = rhs%norm(ChiDG_COMM)


            ! Update smoother(preconditioner), PRIOR to ptc contribution
            !
            ! Mavriplis, "A residual smoothing strategy for accelerating Newton method continuation", 2018.
            !
            if (self%smooth) then
                if (controller%update_preconditioner(data%sdata%lhs,smoother,ChiDG_COMM)) call smoother%update(data%sdata%lhs,data%sdata%rhs)
            end if


            ! Residual-smoothing
            if (self%smooth .and. self%ptc) then
                f_smooth = self%smooth_relax*self%apply_residual_smoother(data,cfln,data%sdata%rhs,smoother,controller)
                b = -rhs-f_smooth
            else if (self%smooth .and. .not. self%ptc) then
                user_msg = 'WARNING: residual smoothing expects pseudo-transient continuation, &
                            but pseudo-transient continuation was found to be turned off. &
                            Residual smoothing is being turned off for consistency.'
                call write_line('------------------------------------------')
                call write_line(user_msg)
                call write_line('------------------------------------------')
            else
                b = -rhs
            end if ! self%smooth


            ! Pseudo-transient continuation
            !   ptc contribution should be before residual smoothing
            !   because the smoother for DG needs the ptc scaling
            !   for stability of the nonlinear smoothing iterations.
            if (self%ptc) then
                call contribute_pseudo_temporal(data,cfln)
            end if ! self%ptc


            ! Convergence check
            call self%record_and_report(resid,timing,niter,cfln(1))
            if (resid < self%tol) exit
            if ( ieee_is_nan(resid) ) then
                user_msg = "newton%solve: NaN residual norm. Check initial solution and operator objects."
                call chidg_signal(FATAL,user_msg)
            end if


            ! Solve system [lhs][dq] = [b] for newton step: [dq]
            call set_forcing_terms(linear_solver)
            call timer_linear%start()
            call linear_solver%solve(lhs,dq,b,preconditioner,controller,data)
            call timer_linear%stop()


            ! Line Search for appropriate step
            q0 = qold
            f0 = resid
            select case (trim(self%search))
                case('Backtrack','backtrack')
                    call self%backtracking(data,system,cfln,q0,qn,f0,fn,fn_fields,f_smooth,niter)
                case('none','')
                    ! Update state, update residual, compute residual norm 
                    qn = q0 + dq
                    data%sdata%q = qn
                    call system%assemble(data,differentiate=NO_DIFF)
                    fn        = rhs%norm(ChiDG_COMM)
                    fn_fields = rhs%norm_fields(ChiDG_COMM)

                case default
                    user_msg = "Invalid nonlinear search algorithm in newton iteration routine. Valid inputs: 'Backtrack', 'none'."
                    call chidg_signal_one(OOPS, user_msg, trim(self%search))
            end select


            ! Accept new solution, clear working storage
            q = qn
            call dq%clear()


            ! Record iteration data
            call self%matrix_iterations%push_back(linear_solver%niter)
            call self%matrix_time%push_back(timer_linear%elapsed())
            call timer_linear%reset()


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

        integer(ik)             :: idom, ielem, itime, ifield
        real(rk)                :: dtau
        real(rk), allocatable   :: mat(:,:)
        type(element_info_t)    :: elem_info

        ! Compute element-local pseudo-timestep
        call compute_pseudo_timestep(data,cfln)

        ! Add mass/dt to sub-block diagonal in dR/dQ
        do idom = 1,data%mesh%ndomains()
            do ielem = 1,data%mesh%domain(idom)%nelem
                elem_info = data%mesh%get_element_info(idom,ielem)
                do itime = 1,data%mesh%domain(idom)%ntime
                    do ifield = 1,data%eqnset(elem_info%eqn_ID)%prop%nprimary_fields()

                        ! get element-local timestep
                        dtau = minval(data%mesh%domain(idom)%elems(ielem)%dtau)
                        mat = data%mesh%domain(idom)%elems(ielem)%mass * (ONE/dtau)
                        call data%sdata%lhs%scale_diagonal(mat,elem_info,ifield,itime)

                    end do !ifield
                end do !itime
            end do !ielem
        end do !idom

        ! Reassemble Matrix
        call data%sdata%lhs%assemble()

    end subroutine contribute_pseudo_temporal
    !*************************************************************************************





    !>  Compute smoothed residual
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   01/15/2019
    !!
    !-------------------------------------------------------------------------------------
    function apply_residual_smoother(self,data,cfln,vector,smoother,controller) result(smoothed_vector)
        class(newton_t),            intent(in)      :: self
        type(chidg_data_t),         intent(inout)   :: data
        real(rk),                   intent(in)      :: cfln(:)
        type(chidg_vector_t),       intent(in)      :: vector
        class(preconditioner_t),    intent(inout)   :: smoother
        class(solver_controller_t), intent(inout)   :: controller

        type(chidg_vector_t)    :: smoothed_vector, scaled_vector, residual, vector_in

        ! Copy input
        vector_in = vector

        ! Apply pseudo-transient scaling: M/dtau
        scaled_vector = pseudo_transient_scaling(data,cfln,vector_in)

        ! Apply smoother: D^{-1}
        smoothed_vector = smoother%apply(data%sdata%lhs,scaled_vector)

    end function apply_residual_smoother
    !************************************************************************************



    !>
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    function pseudo_transient_scaling(data,cfln,vector) result(scaled_vector)
        type(chidg_data_t),     intent(inout)   :: data
        real(rk),               intent(in)      :: cfln(:)
        type(chidg_vector_t),   intent(inout)   :: vector

        integer(ik)             :: idom, ielem, eqn_ID, itime, ifield
        real(rk)                :: dtau
        real(rk),   allocatable :: field(:)
        type(chidg_vector_t)    :: scaled_vector
        type(element_info_t)    :: elem_info

        ! Compute element-local pseudo-timestep
        call compute_pseudo_timestep(data,cfln)

        ! Initialize storage
        scaled_vector = vector

        ! Scale vector by (M/dtau)
        do idom = 1,data%mesh%ndomains()
            do ielem = 1,data%mesh%domain(idom)%nelem

                elem_info = data%mesh%get_element_info(idom,ielem)

                do itime = 1,elem_info%ntime
                    do ifield = 1,data%eqnset(elem_info%eqn_ID)%prop%nprimary_fields()

                        ! get element-local timestep
                        dtau = minval(data%mesh%domain(idom)%elems(ielem)%dtau)

                        ! Retrieve field
                        !field = scaled_vector%get_field(elem_info,ifield,itime)
                        field = vector%get_field(elem_info,ifield,itime)

                        ! Scale field by (M/dtau)
                        field = matmul(data%mesh%domain(idom)%elems(ielem)%mass/dtau,field)

                        ! Store scaled field 
                        call scaled_vector%set_field(field,elem_info,ifield,itime)

                    end do !ifield
                end do !itime
            end do !ielem
        end do !idom

        ! Reassemble
        call scaled_vector%assemble()

    end function pseudo_transient_scaling
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
    subroutine backtracking(self,data,system,cfln,q0,qn,f0,fn,fn_fields,f_smooth,niter)
        class(newton_t),            intent(inout)   :: self
        type(chidg_data_t),         intent(inout)   :: data
        class(system_assembler_t),  intent(inout)   :: system
        real(rk),                   intent(inout)   :: cfln(:)
        type(chidg_vector_t),       intent(inout)   :: q0
        type(chidg_vector_t),       intent(inout)   :: qn
        real(rk),                   intent(inout)   :: f0
        real(rk),                   intent(inout)   :: fn
        real(rk),   allocatable,    intent(inout)   :: fn_fields(:)
        type(chidg_vector_t),       intent(inout)   :: f_smooth
        integer(ik),                intent(inout)   :: niter

        real(rk)                :: alpha, rhs_norm, fn_prev
        real(rk),   allocatable :: rhs_norm_fields(:)
        logical                 :: searching
        integer(ik)             :: step, icfl

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
            call system%assemble(data,differentiate=NO_DIFF)

!           It could be reasonable to track the residual accounting for residual smoothing and ptc as described by 
!           Mavriplis (2018), but it introduces some inconsistencies in other areas so we elect to track the residual of 
!           the original function. Hence these contributions are commented out.
!            if (self%ptc)    data%sdata%rhs = data%sdata%rhs + alpha*pseudo_transient_scaling(data,cfln,data%sdata%dq)
!            if (self%smooth) data%sdata%rhs = data%sdata%rhs + f_smooth

            ! Compute n-th residual norm
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
!            else if (fn > f0) then
!                searching = .true.
            else
                searching = .false.
            end if

            call write_line("       Rn(Q) = ", fn, io_proc=GLOBAL_MASTER, silence=(verbosity<3))

            step = step + 1
        end do

        ! Compute norm by field
        fn_fields = data%sdata%rhs%norm_fields(ChiDG_COMM)


        ! Update cfl
        if (alpha > 0.75_rk) then
            cfln = cfln*self%cfl_up
        else if (alpha < 0.2) then
            cfln = cfln*self%cfl_down
        else
            ! If inbetween, cfl stays the same
        end if


        ! If cfl_max is > 0, enforce
        do icfl = 1,size(cfln)
            if (self%cfl_max > 0. .and. cfln(icfl) > self%cfl_max) cfln(icfl) = self%cfl_max
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

        if ( resid < self%tol ) then
            call write_line('|R| = ', resid, io_proc=GLOBAL_MASTER, silence=(verbosity<4))
            call self%residual_norm%push_back(resid)
            call self%residual_time%push_back(timing)
            call self%matrix_iterations%push_back(0)
            call self%matrix_time%push_back(0._rk)
            call self%newton_iterations%push_back(niter)
            call write_line(niter, resid, cfln, 0, delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))
        ! Compute and store first residual norm 
        else if (niter == 1) then
            call write_line('|R| = ', resid, io_proc=GLOBAL_MASTER, silence=(verbosity<4))
            call self%residual_norm%push_back(resid)
            call self%residual_time%push_back(timing)
            call self%matrix_iterations%push_back(0)
            call self%matrix_time%push_back(0._rk)
            call write_line(niter, resid, cfln, 0, delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))
        end if

    end subroutine record_and_report
    !***********************************************************************************






end module type_newton






