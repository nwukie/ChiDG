module type_newton_subiter
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
    type, extends(nonlinear_solver_t), public :: newton_subiter_t

    contains
        procedure   :: solve
        procedure   :: backtracking
        procedure   :: record_and_report
        procedure   :: update_cfl
!        procedure   :: apply_residual_smoother
    end type newton_subiter_t
    !******************************************************************************************


contains


    !>  Solve for update 'dq'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine solve(self,data,system,linear_solver,preconditioner,solver_controller)
        class(newton_subiter_t),                        intent(inout)           :: self
        type(chidg_data_t),                     intent(inout)           :: data
        class(system_assembler_t),  optional,   intent(inout)           :: system
        class(linear_solver_t),     optional,   intent(inout)           :: linear_solver
        class(preconditioner_t),    optional,   intent(inout),  target  :: preconditioner
        class(solver_controller_t), optional,   intent(inout),  target  :: solver_controller

        character(:),   allocatable :: user_msg
        integer(ik)             :: itime, niter, ierr, icfl
        real(rk)                :: cfl, timing, resid, resid_prev, resid0, resid_new,    &
                                   alpha, f0, fn,fn_pure, forcing_term, residual_ratio
        real(rk), allocatable   :: cfln(:), rnorm0(:), rnorm(:), fn_fields(:)
        type(chidg_vector_t)    :: b, qn, qold, q0, f_smooth, q_sub_iterate, dq_prev, rhs_ptc
        logical                 :: absolute_convergence, relative_convergence, stop_run, iteration_convergence

        type(solver_controller_t),  target  :: default_controller
        class(solver_controller_t), pointer :: controller
        class(preconditioner_t),    pointer :: smoother => null()
        type(precon_jacobi_t),      target  :: jacobi

        integer(ik) :: nsub_iter, max_sub_iter
        logical     :: sub_iteration_convergence, stop_sub_iteration
        
        real(rk)    :: alpha_out, resid_pure
      

        ! Default controller
        controller => default_controller
        if (present(solver_controller)) controller => solver_controller

      
        call write_line('NONLINEAR SOLVER', io_proc=GLOBAL_MASTER, silence=(verbosity<2))
        call write_line("iter","|R(Q)|","CFL", "Linear Solver(niter)", "LHS Updated", "Preconditioner Updated", &
                        delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))


        ! start timer
        call self%timer%reset()
        call self%timer%start()


        ! Initialize smoother
        select case (self%smoother)
            case('preconditioner','Preconditioner')
                smoother => preconditioner
            case('jacobi','Jacobi')
                smoother => jacobi
                call smoother%init(data)
            case default
                call chidg_signal(FATAL,"newton_subiter%solve: invalid smoother. 'default' or 'jacobi'.")
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

        max_sub_iter = 3 


        ! Nonlinear iteration, while not converged
        do while ( absolute_convergence .and. &
                   relative_convergence .and. &
                   iteration_convergence .and. &
                   (.not. stop_run) )
            niter = niter + 1

            ! Store the value of the current inner iteration solution (k) 
            ! for the solution update (n+1), q_(n+1)_k
            qold = q

            dq_prev = ZERO*dq

           nsub_iter=0
           sub_iteration_convergence = .true.
           stop_sub_iteration = .false.

            do while (sub_iteration_convergence .and. &
                    (.not. stop_sub_iteration))

            !do while ((.not. stop_sub_iteration))
                nsub_iter = nsub_iter + 1
                call write_line('SUBITERATION ',nsub_iter, io_proc=GLOBAL_MASTER, silence=(verbosity<2))

                q_sub_iterate = q
                dq_prev = q_sub_iterate - qold


                ! Update Spatial Residual and Linearization (rhs, lin)
                residual_ratio = resid/resid_prev
                if ( niter <= 2) residual_ratio = ONE
                call write_line('ASSEMBLE...', io_proc=GLOBAL_MASTER, silence=(verbosity<2))
                call system%assemble( data,             &
                                      timing=timing,    &
                                      differentiate=controller%update_lhs(lhs,niter,residual_ratio) )



                if ((niter == 1) .and. (nsub_iter==1)) then
                    rnorm0 = rhs%norm_fields(ChiDG_COMM)
                    cfln = rnorm0
                    cfln = self%cfl0
                end if

                resid_pure = rhs%norm(ChiDG_COMM)


                ! Pseudo-transient continuation
                !   ptc contribution should be before residual smoothing
                !   because the smoother for DG needs the ptc scaling
                !   for stability of the nonlinear smoothing iterations.
                if (self%ptc) then
                    call contribute_pseudo_temporal(data,cfln)
                    rhs_ptc = pseudo_transient_scaling(data,cfln,dq_prev)
                end if ! self%ptc

                if ((niter == 1) .and. (nsub_iter==1)) then
                    resid0 = rhs%norm(ChiDG_COMM)
                    rnorm0 = rhs%norm_fields(ChiDG_COMM)
                end if

                ! Residual-smoothing
!                if (self%smooth .and. self%ptc) then
!                    f_smooth = self%smooth_relax*self%apply_residual_smoother(data,cfln,data%sdata%rhs,smoother,controller)
!                    b = -rhs-f_smooth
!                else if (self%smooth .and. .not. self%ptc) then
!                    user_msg = 'WARNING: residual smoothing expects pseudo-transient continuation, &
!                                but pseudo-transient continuation was found to be turned off. &
!                                Residual smoothing is being turned off for consistency.'
!                    call write_line('------------------------------------------')
!                    call write_line(user_msg)
!                    call write_line('------------------------------------------')
!                else
!                    b = -rhs - rhs_ptc
!                end if ! self%smooth

                if (self%ptc) then
                    b = -rhs - rhs_ptc
                else
                    b = -rhs
                end if



                ! Compute and store residual norm for each field
                resid_prev = resid
                resid      = b%norm(ChiDG_COMM)
                rnorm = b%norm_fields(ChiDG_COMM)


                ! During the first few iterations, allow the initial residual norm to update
                ! if it has increased. Otherwise, if a solution converged to an error floor
                ! and the residual raised a little bit, then the CFL would essentially reset
                ! from infinity to CFL0, which we do not want.
                if (niter < 5) then
                    where (rnorm > rnorm0)
                        rnorm0 = rnorm0 + 0.3_rk*(rnorm - rnorm0)
                    end where
                end if


                ! Convergence check
                call self%record_and_report(resid_pure,timing,niter, nsub_iter,cfln(1))
                if (resid < self%tol) exit



                if ( ieee_is_nan(resid) ) then
                    user_msg = "newton_subiter%solve: NaN residual norm. Check initial solution and operator objects."
                    call chidg_signal(FATAL,user_msg)
                end if




                ! Solve system [lhs][dq] = [b] for newton_subiter step: [dq]
                call set_forcing_terms(linear_solver)
                call write_line('LINEAR SOLVER...', io_proc=GLOBAL_MASTER, silence=(verbosity<2))
                call linear_solver%solve(lhs,dq,b,preconditioner,controller,data)


                ! Line Search for appropriate step
                q0 = qold
                f0 = resid
                call write_line('BACKTRACKING...', io_proc=GLOBAL_MASTER, silence=(verbosity<2))
                select case (trim(self%search))
                    case('Backtrack','backtrack')
                        call self%backtracking(data,system,cfln,q0,qn,q_sub_iterate,f0,fn,fn_pure,fn_fields,f_smooth, alpha_out)
                    case('none','')
                        ! Update state, update residual, compute residual norm 
                        qn = q_sub_iterate+ dq
                        data%sdata%q = qn
                        call system%assemble(data,differentiate=NO_DIFF)

                        if (self%ptc)    rhs = rhs + pseudo_transient_scaling(data,cfln,data%sdata%dq)
                        if (self%smooth) rhs = rhs + f_smooth
                        fn        = rhs%norm(ChiDG_COMM)
                        fn_fields = rhs%norm_fields(ChiDG_COMM)

                    case default
                        user_msg = "Invalid nonlinear search algorithm in newton iteration routine. Valid inputs: 'Backtrack', 'none'."
                        call chidg_signal_one(OOPS, user_msg, trim(self%search))
                end select


                ! Accept new solution, update cfl using new residual, clear working storage
                q = qn
                !dq_prev = dq
                call dq%clear()

                !stop_sub_iteration = (nsub_iter < max_sub_iter)
                stop_sub_iteration = (max_sub_iter <= nsub_iter)

                !sub_iteration_convergence  = (fn/resid0 > self%rtol)

                sub_iteration_convergence = (fn > 1.0)


                call write_line(niter+1, fn_pure, cfln(1), linear_solver%niter, controller%lhs_updated, controller%preconditioner_updated, delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))
            end do !nsub_iter

            ! Update cfl
            if (alpha_out > 0.75_rk) then
                cfln = cfln*self%cfl_up
            else if (alpha_out< 0.1) then
                cfln = cfln*self%cfl_down
            else
                ! If inbetween, cfl stays the same
            end if


            ! If cfl_max is > 0, enforce
            do icfl = 1,size(cfln)
                if (self%cfl_max > 0. .and. cfln(icfl) > self%cfl_max) cfln(icfl) = self%cfl_max
            end do


            ! Record iteration data
            !call self%update_cfl(rnorm0,fn_fields,cfln)
            call self%matrix_iterations%push_back(linear_solver%niter)
            call self%matrix_time%push_back(linear_solver%timer%elapsed())
            call self%residual_norm%push_back(fn_pure)
            call self%residual_time%push_back(timing)
            absolute_convergence  = (fn_pure > self%tol)
            relative_convergence  = (fn_pure/resid0 > self%rtol)
            iteration_convergence = (niter < self%nmax) .or. self%nmax <= 0
            inquire(file='STOP', exist=stop_run)



            ! Print iteration information
            call write_line(niter+1, fn_pure, cfln(1), linear_solver%niter, controller%lhs_updated, controller%preconditioner_updated, delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))
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
        class(newton_subiter_t),            intent(in)      :: self
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
                cfln = self%cfl0*(rnorm0(1)/rnorm(1))
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

        ! If cfl_max is > 0, enforce
        do icfl = 1,size(cfln)
            if (self%cfl_max > 0. .and. cfln(icfl) > self%cfl_max) cfln(icfl) = self%cfl_max
        end do

        ! Strategy for cfl across fields
        select case(trim(self%cfl_fields))
            case('minimum')
                cfln = minval(cfln)
            case('average')
                cfln = sum(cfln)/size(cfln)
            case default
                call chidg_signal(FATAL,"newton_subiter%update_cfl: invalid selection of cfl_fields behavior: 'minimum' or 'average'.")
        end select


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



















!    !>  Compute smoothed residual
!    !!
!    !!  @author Nathan A. Wukie(AFRL)
!    !!  @date   01/15/2019
!    !!
!    !-------------------------------------------------------------------------------------
!    function apply_residual_smoother(self,data,cfln,vector,smoother,controller) result(scaled_vector)
!        class(newton_subiter_t),            intent(in)      :: self
!        type(chidg_data_t),         intent(inout)   :: data
!        real(rk),                   intent(in)      :: cfln(:)
!        type(chidg_vector_t),       intent(in)      :: vector
!        class(preconditioner_t),    intent(inout)   :: smoother
!        class(solver_controller_t), intent(inout)   :: controller
!
!        type(chidg_vector_t)    :: smoothed_vector, scaled_vector, residual
!        integer(ik) :: ismooth
!
!        ! Update smoother(preconditioner)
!        if (controller%update_preconditioner(data%sdata%lhs,smoother)) call smoother%update(data%sdata%lhs,vector)
!
!        ! Apply residual smoothing
!        smoothed_vector = ZERO*vector
!        do ismooth = 1,self%nsmooth
!            residual = vector - chidg_mv(data%sdata%lhs,smoothed_vector)
!            smoothed_vector = smoother%apply(data%sdata%lhs,residual)
!        end do
!
!        ! Apply pseudo-transient scaling
!        scaled_vector = pseudo_transient_scaling(data,cfln,smoothed_vector)
!
!    end function apply_residual_smoother
!    !************************************************************************************


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
    subroutine backtracking(self,data,system,cfln,q0,qn,q_sub,f0,fn,fn_pure,fn_fields,f_smooth, alpha_out)
        class(newton_subiter_t),            intent(inout)   :: self
        type(chidg_data_t),         intent(inout)   :: data
        class(system_assembler_t),  intent(inout)   :: system
        real(rk),                   intent(inout)   :: cfln(:)
        type(chidg_vector_t),       intent(inout)   :: q0
        type(chidg_vector_t),       intent(inout)   :: qn
        type(chidg_vector_t),       intent(inout)   :: q_sub
        real(rk),                   intent(inout)   :: f0
        real(rk),                   intent(inout)   :: fn
        real(rk),                   intent(inout)   :: fn_pure
        real(rk),   allocatable,    intent(inout)   :: fn_fields(:)
        type(chidg_vector_t),       intent(inout)   :: f_smooth
        real(rk),                   intent(inout)   :: alpha_out

        real(rk)                :: alpha, rhs_norm, fn_prev
        real(rk),   allocatable :: rhs_norm_fields(:)
        logical                 :: searching
        integer(ik)             :: step, icfl
        type(chidg_vector_t)    :: dq_sub

        step = 0
        searching = .true.
        call write_line("Initial residual norm: ", f0, io_proc=GLOBAL_MASTER, silence=(verbosity<2))
        do while (searching)


            call write_line("BACKTRACKING STEP ", step+1, io_proc=GLOBAL_MASTER, silence=(verbosity<2))
            alpha = TWO**(-real(step,rk)) 
            call write_line("       Testing newton_subiter direction with 'alpha' = ", alpha, io_proc=GLOBAL_MASTER, silence=(verbosity<3))

            ! Advance solution along newton direction
            qn = q_sub + alpha*data%sdata%dq

            ! Set working solution. Test residual at (q). Do not differentiate
            data%sdata%q = qn

            ! Compute new function value and norm
            call system%assemble(data,differentiate=NO_DIFF)

            fn_pure = data%sdata%rhs%norm(ChiDG_COMM)
            ! Add globalization contributions
            dq_sub = qn-q0
            if (self%ptc)    data%sdata%rhs = data%sdata%rhs + pseudo_transient_scaling(data,cfln,dq_sub)
            !if (self%smooth) data%sdata%rhs = data%sdata%rhs + f_smooth

            ! Compute n-th residual norm
            fn = data%sdata%rhs%norm(ChiDG_COMM)
            call write_line("Tested residual norm: ", fn, io_proc=GLOBAL_MASTER, silence=(verbosity<2))

            ! Test for |R| increasing too much or NaN. 
            !   If residual is reasonably large, still allow some growth.
            !   If the residual is small enough, we don't want any growth.
            if (ieee_is_nan(fn)) then
                searching = .true.
            !else if ( (fn > 1.e-3_rk) .and. (fn > 2.0_rk*f0) ) then
            !    searching = .true.
            !else if ( (fn < 1.e-3_rk) .and. (fn > f0) ) then
            !    searching = .true.
            else if (fn > f0) then
                searching = .true.
            else
                searching = .false.
                call write_line("     Accepted  'alpha' = ", alpha, io_proc=GLOBAL_MASTER, silence=(verbosity<2))
                call write_line("     |PTC-R|=", fn, io_proc=GLOBAL_MASTER, silence=(verbosity<2))
            end if

            call write_line("       Rn(Q) = ", fn, io_proc=GLOBAL_MASTER, silence=(verbosity<3))

            step = step + 1
        end do

        ! Compute norm by field
        fn_fields = data%sdata%rhs%norm_fields(ChiDG_COMM)

        alpha_out = alpha

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
    subroutine record_and_report(self,resid,timing,niter,nsub_iter,cfln)
        class(newton_subiter_t),    intent(inout)   :: self
        real(rk),           intent(in)      :: resid
        real(rk),           intent(in)      :: timing
        integer(ik),        intent(in)      :: niter
        integer(ik),        intent(in)      :: nsub_iter
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
        else if ((niter == 1) .and. (nsub_iter==1)) then
            call write_line('|R| = ', resid, io_proc=GLOBAL_MASTER, silence=(verbosity<4))
            call self%residual_norm%push_back(resid)
            call self%residual_time%push_back(timing)
            call self%matrix_iterations%push_back(0)
            call self%matrix_time%push_back(0._rk)
            call write_line(niter, resid, cfln, 0, delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))
        end if

    end subroutine record_and_report
    !***********************************************************************************






end module type_newton_subiter






