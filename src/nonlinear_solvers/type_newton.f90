module type_newton
#include <messenger.h>
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"
#include "petsc/finclude/petscmat.h"
    use petscmat,               only: MatSetValues, ADD_VALUES



    use petscksp,               only: tKSP, KSPCreate, KSPSolve, KSPSetOperators, KSPSetType, KSPGetPC, KSPSetUp, &
                                      KSPSetTolerances, PETSC_DEFAULT_INTEGER, PETSC_DEFAULT_REAL, KSPGMRESSetRestart, &
                                      KSPGMRESSetCGSRefinementType, KSP_GMRES_CGS_REFINE_ALWAYS, KSP_GMRES_CGS_REFINE_IFNEEDED, &
                                      KSPGetIterationNumber

    use petscpc,                only: tPC, PCSetType, PCFactorSetLevels, PCFactorSetShiftType, MAT_SHIFT_POSITIVE_DEFINITE
!    use petscpc,                only: PCSetType, tPC, PCFactorSetLevels, PCGAMGAGG, PCHYPRE
!    use petscpc,                only: PCHYPRE, PCSetType

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
    use type_timer,             only: timer_t
    use type_element_info,      only: element_info_t
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
        procedure   :: update_cfl
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
        integer(ik)             :: itime, niter, ierr, icfl
        real(rk)                :: cfl, timing, resid, resid_prev, resid0, resid_new,    &
                                   alpha, f0, fn, forcing_term, residual_ratio, testval
        real(rk), allocatable   :: cfln(:), rnorm0(:), rnorm(:), fn_fields(:)
        type(chidg_vector_t)    :: b, qn, qold, q0, f_smooth, test
        logical                 :: absolute_convergence, relative_convergence, stop_run, iteration_convergence

        type(solver_controller_t),  target  :: default_controller
        class(solver_controller_t), pointer :: controller
        class(preconditioner_t),    pointer :: smoother => null()
        type(precon_jacobi_t),      target  :: jacobi
        type(timer_t)                       :: timer_linear

        PetscErrorCode :: perr

        KSP :: ksp
        PC  :: pc


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


                call KSPCreate(ChiDG_COMM%mpi_val,ksp,perr)
                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling KSPCreate.')

                call KSPSetOperators(ksp,data%sdata%lhs%petsc_matrix,data%sdata%lhs%petsc_matrix,perr)
                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling KSPSetOperators.')

                call KSPSetType(ksp,KSPGMRES,perr)
                !call KSPSetType(ksp,KSPDGMRES,perr)
                !call KSPSetType(ksp,KSPBCGS,perr)
                !call KSPSetType(ksp,KSPGCR,perr)
                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling KSPSetType.')

                
                !*******    Preconditioners   *********!
                call KSPGetPC(ksp,pc,perr)
                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling KSPGetPC.')

!                call PCSetType(pc,'none',perr)
!                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling PCSetType.')

!                call PCSetType(pc,'gamg',perr)
!                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling PCSetType.')
!                call PCGAMGSetType(gc, PCGAMGAGG, perr)
!                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling PCGAMGSetType.')

!                call PCSetType(pc, PCHYPRE, perr)
!                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling PCSetType.')
!                call PCHYPRESetType(pc,'boomeramg', perr)
!                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling PCHYPRESetType.')
!                call PCHYPRESetType(pc,'parasails', perr)
!                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling PCHYPRESetType.')

!                call PCSetType(pc, PCSOR, perr)
!                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling PCSetType.')

!                call PCSetType(pc,'jacobi',perr)
!                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling PCSetType.')

!                call PCSetType(pc,'ilu',perr)
!                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling PCSetType.')
!                call PCFactorSetShiftType(pc, MAT_SHIFT_POSITIVE_DEFINITE, perr)
!                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling PCFactorSetShiftType.')
!                call PCFactorSetLevels(pc,0,perr)
!                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling PCFactorSetLevels.')

!                call PCSetType(pc,PCSPAI,perr)
!                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling PCSetType.')
                !*******    Preconditioners   *********!






                call KSPSetFromOptions(ksp,perr)
                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling KSPSetFromOptions.')

                call KSPSetTolerances(ksp, linear_solver%rtol, linear_solver%tol, PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, perr)
                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling KSPSetTolerances.')

                call KSPGMRESSetRestart(ksp, linear_solver%nkrylov, perr)
                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling KSPGMRESSetRestart.')
                
                call KSPGMRESSetCGSRefinementType(ksp, KSP_GMRES_CGS_REFINE_ALWAYS, perr)
                !call KSPGMRESSetCGSRefinementType(ksp, KSP_GMRES_CGS_REFINE_IFNEEDED, perr)
                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling KSPGMRESCGSSetRefinementType.')

                call KSPSetUp(ksp,perr)
                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling KSPSetUp.')

            else

                call KSPSetOperators(ksp,data%sdata%lhs%petsc_matrix,data%sdata%lhs%petsc_matrix,perr)
                if (perr /= 0) call chidg_signal(FATAL,'newton: error calling KSPSetOperators.')

            end if
























            if (niter == 1) then
                rnorm0 = rhs%norm_fields(ChiDG_COMM)
                cfln = rnorm0
                cfln = self%cfl0
            end if



            ! Pseudo-transient continuation
            !   ptc contribution should be before residual smoothing
            !   because the smoother for DG needs the ptc scaling
            !   for stability of the nonlinear smoothing iterations.
            if (self%ptc) then
                call contribute_pseudo_temporal(data,cfln)
            end if ! self%ptc


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



            ! Compute and store residual norm for each field
            resid_prev = resid
            resid      = b%norm(ChiDG_COMM)
            if (niter == 1) then
                resid0 = b%norm(ChiDG_COMM)
                rnorm0 = b%norm_fields(ChiDG_COMM)
            end if
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
            call self%record_and_report(resid,timing,niter,cfln(1))
            if (resid < self%tol) exit



            if ( ieee_is_nan(resid) ) then
                user_msg = "newton%solve: NaN residual norm. Check initial solution and operator objects."
                call chidg_signal(FATAL,user_msg)
            end if


            ! Solve system [lhs][dq] = [b] for newton step: [dq]
            call set_forcing_terms(linear_solver)

























            !call linear_solver%solve(lhs,dq,b,preconditioner,controller,data)

            call timer_linear%start()
            call KSPSolve(ksp,b%petsc_vector,dq%petsc_vector,perr)
            call timer_linear%stop()
            if (perr /= 0) call chidg_signal(FATAL,'newton: error calling KSPSolve.')
            call KSPGetIterationNumber(ksp, linear_solver%niter, perr)
            if (perr /= 0) call chidg_signal(FATAL,'newton: error calling KSPGetIterationNumber.')



            ! Line Search for appropriate step
            q0 = qold
            f0 = resid
            select case (trim(self%search))
                case('Backtrack','backtrack')
                    call self%backtracking(data,system,cfln,q0,qn,f0,fn,fn_fields,f_smooth)
                case('none','')
                    ! Update state, update residual, compute residual norm 
                    qn = q0 + dq
                    data%sdata%q = qn
                    call system%assemble(data,differentiate=.false.)

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
            !call self%update_cfl(rnorm0,fn_fields,cfln)
            call dq%clear()


            ! Record iteration data
            call self%matrix_iterations%push_back(linear_solver%niter)
            !call self%matrix_time%push_back(linear_solver%timer%elapsed())
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
                call chidg_signal(FATAL,"newton%update_cfl: invalid selection of cfl_fields behavior: 'minimum' or 'average'.")
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

        integer(ik)                 :: idom, ielem, eqn_ID, itime, ifield, nterms, rstart, rend, cstart, cend, imat, dof_start, nrows, ncols, iarray, row_index_start, col_index_start, i
        integer(ik), allocatable, dimension(:) :: col_indices, row_indices
        real(rk)                    :: dtau
        real(rk), allocatable       :: mat(:,:)

        PetscErrorCode :: ierr


        ! Compute element-local pseudo-timestep
        call compute_pseudo_timestep(data,cfln)


        ! Add mass/dt to sub-block diagonal in dR/dQ
        if (data%sdata%q%petsc_vector_created) then

            do idom = 1,data%mesh%ndomains()
                do ielem = 1,data%mesh%domain(idom)%nelem
                    dof_start = data%mesh%domain(idom)%elems(ielem)%dof_start
                    eqn_ID = data%mesh%domain(idom)%elems(ielem)%eqn_ID
                    do itime = 1,data%mesh%domain(idom)%ntime
                        do ifield = 1,data%eqnset(eqn_ID)%prop%nprimary_fields()

                            ! get element-local timestep
                            if (size(cfln) == 1) then
                                dtau = data%mesh%domain(idom)%elems(ielem)%dtau(1)
                            else
                                dtau = data%mesh%domain(idom)%elems(ielem)%dtau(ifield)
                            end if

                            ! Need to compute row and column extends in diagonal so we can
                            ! selectively apply the mass matrix to the sub-block diagonal
                            nterms = data%mesh%domain(idom)%elems(ielem)%nterms_s


                            row_index_start = dof_start + (ifield-1)*nterms + (data%eqnset(eqn_ID)%prop%nprimary_fields()*nterms)*(itime-1)
                            row_indices     = [(i, i=row_index_start,(row_index_start+nterms-1),1)]

                            ! Diagonal addition so column indices are identical to row indices
                            col_index_start = row_index_start
                            col_indices     = row_indices


                            mat = data%mesh%domain(idom)%elems(ielem)%mass * (ONE/dtau)

                            nrows = 1
                            ncols = size(mat,2)
                            do iarray = 1,size(mat,1)
                                ! subtract 1 from indices since petsc is 0-based
                                call MatSetValues(data%sdata%lhs%petsc_matrix,nrows,[row_index_start + (iarray-1) - 1],ncols,col_indices-1,mat(iarray,:),ADD_VALUES,ierr)
                                if (ierr /= 0) call chidg_signal(FATAL,"chidg_matrix%petsc_store: error calling MatSetValues.")
                            end do 


                        end do !ifield
                    end do !itime
                end do !ielem
            end do !idom

            ! Reassemble Matrix
            call data%sdata%lhs%assemble()

            ! Update stamp
            call date_and_time(values=data%sdata%lhs%stamp)


        else

            do idom = 1,data%mesh%ndomains()
                do ielem = 1,data%mesh%domain(idom)%nelem
                    eqn_ID = data%mesh%domain(idom)%elems(ielem)%eqn_ID
                    do itime = 1,data%mesh%domain(idom)%ntime
                        do ifield = 1,data%eqnset(eqn_ID)%prop%nprimary_fields()

                            ! get element-local timestep
                            dtau = data%mesh%domain(idom)%elems(ielem)%dtau(ifield)

                            ! Need to compute row and column extends in diagonal so we can
                            ! selectively apply the mass matrix to the sub-block diagonal
                            nterms = data%mesh%domain(idom)%elems(ielem)%nterms_s
                            rstart = 1 + (ifield-1) * nterms
                            rend   = (rstart-1) + nterms
                            cstart = rstart                 ! since it is square
                            cend   = rend                   ! since it is square

                            ! Add mass matrix divided by dt to the block diagonal
                            imat = data%sdata%lhs%dom(idom)%lblks(ielem,itime)%get_diagonal()
                            data%sdata%lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(rstart:rend,cstart:cend)  =  data%sdata%lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(rstart:rend,cstart:cend)  +  data%mesh%domain(idom)%elems(ielem)%mass*(ONE/dtau)


                        end do !ifield
                    end do !itime
                end do !ielem
            end do !idom

            ! Update stamp
            call date_and_time(values=data%sdata%lhs%stamp)

        end if

    end subroutine contribute_pseudo_temporal
    !*************************************************************************************







    !>  Compute smoothed residual
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   01/15/2019
    !!
    !-------------------------------------------------------------------------------------
    function apply_residual_smoother(self,data,cfln,vector,smoother,controller) result(scaled_vector)
        class(newton_t),            intent(in)      :: self
        type(chidg_data_t),         intent(inout)   :: data
        real(rk),                   intent(in)      :: cfln(:)
        type(chidg_vector_t),       intent(in)      :: vector
        class(preconditioner_t),    intent(inout)   :: smoother
        class(solver_controller_t), intent(inout)   :: controller

        type(chidg_vector_t)    :: smoothed_vector, scaled_vector, residual
        integer(ik) :: ismooth

        ! Update smoother(preconditioner)
        if (controller%update_preconditioner(data%sdata%lhs,smoother)) call smoother%update(data%sdata%lhs,vector)

        ! Apply residual smoothing
        smoothed_vector = ZERO*vector
        do ismooth = 1,self%nsmooth
            residual = vector - chidg_mv(data%sdata%lhs,smoothed_vector)
            smoothed_vector = smoother%apply(data%sdata%lhs,residual)
        end do

        ! Apply pseudo-transient scaling
        scaled_vector = pseudo_transient_scaling(data,cfln,smoothed_vector)

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
        type(element_info_t)    :: element_info

        ! Compute element-local pseudo-timestep
        call compute_pseudo_timestep(data,cfln)

        ! Initialize storage
        scaled_vector = vector

        ! Scale vector by (M/dtau)
        do idom = 1,data%mesh%ndomains()
            do ielem = 1,data%mesh%domain(idom)%nelem
                eqn_ID = data%mesh%domain(idom)%elems(ielem)%eqn_ID

                element_info = element_info_t(idomain_g  = data%mesh%domain(idom)%elems(ielem)%idomain_g,    &
                                              idomain_l  = data%mesh%domain(idom)%elems(ielem)%idomain_l,    &
                                              ielement_g = data%mesh%domain(idom)%elems(ielem)%ielement_g,   &
                                              ielement_l = data%mesh%domain(idom)%elems(ielem)%ielement_l,   &
                                              iproc      = data%mesh%domain(idom)%elems(ielem)%iproc,        &
                                              pelem_ID   = NO_ID,                                            &
                                              eqn_ID     = data%mesh%domain(idom)%elems(ielem)%eqn_ID,       &
                                              nfields    = data%mesh%domain(idom)%elems(ielem)%neqns,        &
                                              nterms_s   = data%mesh%domain(idom)%elems(ielem)%nterms_s,     &
                                              nterms_c   = data%mesh%domain(idom)%elems(ielem)%nterms_c,     &
                                              dof_start  = data%mesh%domain(idom)%elems(ielem)%dof_start)


                do itime = 1,data%mesh%domain(idom)%ntime
                    do ifield = 1,data%eqnset(eqn_ID)%prop%nprimary_fields()

                        ! get element-local timestep
                        if (size(cfln) == 1) then
                            dtau = data%mesh%domain(idom)%elems(ielem)%dtau(1)
                        else
                            dtau = data%mesh%domain(idom)%elems(ielem)%dtau(ifield)
                        end if


                        ! Retrieve field
                        field = scaled_vector%get_field(element_info,ifield,itime)

                        ! Scale field by (M/dtau)
                        field = matmul(data%mesh%domain(idom)%elems(ielem)%mass/dtau,field)

                        ! Store scaled field 
                        call scaled_vector%set_field(field,element_info,ifield,itime)

                    end do !ifield
                end do !itime
            end do !ielem
        end do !idom

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
    subroutine backtracking(self,data,system,cfln,q0,qn,f0,fn,fn_fields,f_smooth)
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
            call system%assemble(data,differentiate=.false.)

            ! Add globalization contributions
!            if (self%ptc)    data%sdata%rhs = data%sdata%rhs + alpha*pseudo_transient_scaling(data,cfln,data%sdata%dq)
            if (self%smooth) data%sdata%rhs = data%sdata%rhs + f_smooth

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
            !else if (fn > f0) then
            !    searching = .true.
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
        else if (alpha < 0.1) then
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






