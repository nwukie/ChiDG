module type_petsc_nonlinear
#include <messenger.h>
#include "petsc/finclude/petscsnes.h"
    use petscsnes,               only: tSNES, SNESCreate

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
    use type_chidg_vector
    use ieee_arithmetic,        only: ieee_is_nan
    implicit none
    private



    !>  Petsc Scalable Nonlinear Equation Solver: SNES
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2019
    !!
    !------------------------------------------------------------------------------------------
    type, extends(nonlinear_solver_t), public :: petsc_nonlinear_t

        SNES    :: snes
        logical :: petsc_initialized = .false.

    contains
        procedure   :: solve
        procedure   :: record_and_report
        procedure   :: initialize_petsc
    end type petsc_nonlinear_t
    !******************************************************************************************


contains


    !>  Solve for update 'dq'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine solve(self,data,system,linear_solver,preconditioner,solver_controller)
        class(petsc_nonlinear_t),               intent(inout)           :: self
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
        type(timer_t)                       :: timer_linear

        PetscErrorCode :: perr

        ! NOT FINISHED!!!
        call chidg_signal(FATAL,'petsc_nonlinear%solve: implementation is incomplete and not currently functional.')

        ! Initialize solver environment
        if (.not. self%petsc_initialized) call self%initialize_petsc()


        !call SNESSetFunction(self%snes, data%sdata%rhs%petsc_vector,



        !call SNESSolve(self%snes, PETSC_NULL_OBJECT, data%sdata%q%petsc_vector, perr)
        !if (perr /= 0) call chidg_signal(FATAL,'petsc_nonlinear%solve: error calling SNESSolve.')





    end subroutine solve
    !*****************************************************************************************




    !>
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine initialize_petsc(self)
        class(petsc_nonlinear_t),   intent(inout)   :: self

        PetscErrorCode :: perr

        call SNESCreate(ChiDG_COMM%mpi_val, self%snes, perr)
        if (perr /= 0) call chidg_signal(FATAL,'petsc_nonlinear: error calling SNESCreate.')


    end subroutine initialize_petsc
    !***********************************************************************************



    !>
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine record_and_report(self,resid,timing,niter,cfln)
        class(petsc_nonlinear_t),   intent(inout)   :: self
        real(rk),                   intent(in)      :: resid
        real(rk),                   intent(in)      :: timing
        integer(ik),                intent(in)      :: niter
        real(rk),                   intent(in)      :: cfln

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






end module type_petsc_nonlinear
