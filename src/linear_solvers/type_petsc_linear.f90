module type_petsc_linear
#include <messenger.h>
#include "petsc/finclude/petscksp.h"
    use petscksp,               only: tKSP, KSPCreate, KSPSolve, KSPSetOperators, KSPSetType, KSPGetPC, KSPSetUp, &
                                      KSPSetTolerances, PETSC_DEFAULT_INTEGER, PETSC_DEFAULT_REAL, KSPGMRESSetRestart, &
                                      KSPGMRESSetCGSRefinementType, KSP_GMRES_CGS_REFINE_ALWAYS, KSP_GMRES_CGS_REFINE_IFNEEDED, &
                                      KSPGetIterationNumber, KSPDestroy

    use mod_kinds,              only: rk, ik
    use mod_chidg_mpi,          only: ChiDG_COMM, GLOBAL_MASTER
    use mod_io,                 only: verbosity
    use mpi_f08

    use type_timer,             only: timer_t
    use type_linear_solver,     only: linear_solver_t 
    use type_preconditioner,    only: preconditioner_t
    use type_solver_controller, only: solver_controller_t
    use type_chidg_data,        only: chidg_data_t
    use type_chidg_vector
    use type_chidg_matrix
    implicit none


    !>  Flexible variant of the Generalized Minimum Residual algorithm for solving nonsymmetric
    !!  systems of linear equations.
    !!
    !!  Algorithms:
    !!      : Number of Krylov vectors       (&linear_solve  nkrylov=2000 /)
    !!      : Orthogonalization              (&linear_solve  orthogonalization='CGS' or 'MGS' /)
    !!      : Inner fgmres iteration   (&linear_solve  inner_fgmres=.true. /)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/11/2016
    !!
    !---------------------------------------------------------------------------------------------
    type, public, extends(linear_solver_t) :: petsc_linear_t

        KSP     :: ksp
        logical :: petsc_initialized = .false.

    contains
        procedure   :: solve
    end type petsc_linear_t
    !*********************************************************************************************


contains





    !> Solution routine
    !!
    !!  NOTE: this subroutine is recursive so local variables should NOT be declared with the
    !!        SAVE attribute. The declaration and allocation of krylov vectors was relocated
    !!        to the object so that their allocation can be stored without using the SAVE 
    !!        attribute.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/11/2016
    !!
    !!  @author Nathan A. Wukie (AFRL) 
    !!  @date   6/23/2016
    !!  @note   parallelization
    !!
    !---------------------------------------------------------------------------------------------
    recursive subroutine solve(self,A,x,b,M,solver_controller,data)
        class(petsc_linear_t),      intent(inout)               :: self
        type(chidg_matrix_t),       intent(inout)               :: A
        type(chidg_vector_t),       intent(inout)               :: x
        type(chidg_vector_t),       intent(inout)               :: b
        class(preconditioner_t),    intent(inout), optional     :: M
        class(solver_controller_t), intent(inout), optional     :: solver_controller
        type(chidg_data_t),         intent(in),    optional     :: data

        PetscErrorCode :: perr

        if (self%petsc_initialized) then
            call KSPSetOperators(self%ksp,A%petsc_matrix,A%petsc_matrix,perr)
            if (perr /= 0) call chidg_signal(FATAL,'petsc_linear: error calling KSPSetOperators.')

        else

            call KSPCreate(ChiDG_COMM%mpi_val,self%ksp,perr)
            if (perr /= 0) call chidg_signal(FATAL,'petsc_linear: error calling KSPCreate.')

            call KSPSetOperators(self%ksp,A%petsc_matrix,A%petsc_matrix,perr)
            if (perr /= 0) call chidg_signal(FATAL,'petsc_linear: error calling KSPSetOperators.')

            call KSPSetType(self%ksp,KSPGMRES,perr)
            if (perr /= 0) call chidg_signal(FATAL,'petsc_linear: error calling KSPSetType.')

            call KSPSetFromOptions(self%ksp,perr)
            if (perr /= 0) call chidg_signal(FATAL,'petsc_linear: error calling KSPSetFromOptions.')

            call KSPSetTolerances(self%ksp, self%rtol, self%tol, PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, perr)
            if (perr /= 0) call chidg_signal(FATAL,'petsc_linear: error calling KSPSetTolerances.')

            call KSPGMRESSetRestart(self%ksp, self%nkrylov, perr)
            if (perr /= 0) call chidg_signal(FATAL,'petsc_linear: error calling KSPGMRESSetRestart.')
            
            call KSPGMRESSetCGSRefinementType(self%ksp, KSP_GMRES_CGS_REFINE_IFNEEDED, perr)
            if (perr /= 0) call chidg_signal(FATAL,'petsc_linear: error calling KSPGMRESCGSSetRefinementType.')

            call KSPSetUp(self%ksp,perr)
            if (perr /= 0) call chidg_signal(FATAL,'petsc_linear: error calling KSPSetUp.')

            self%petsc_initialized = .true.

        end if



        ! Solve
        call KSPSolve(self%ksp,b%wrapped_petsc_vector%petsc_vector,x%wrapped_petsc_vector%petsc_vector,perr)
        if (perr /= 0) call chidg_signal(FATAL,'petsc_linear: error calling KSPSolve.')


        ! Get info
        call KSPGetIterationNumber(self%ksp, self%niter, perr)
        if (perr /= 0) call chidg_signal(FATAL,'petsc_linear: error calling KSPGetIterationNumber.')


    end subroutine solve
    !************************************************************************************************************



    !>  Tear down activities.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2019
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine tear_down(self)
        class(petsc_linear_t),  intent(inout)   :: self

        PetscErrorCode :: perr

        call KSPDestroy(self%ksp,perr)
        if (perr /= 0) call chidg_signal(FATAL,'petsc_linear%tear_down: error calling KSPDestroy.')
        self%petsc_initialized = .false.

    end subroutine tear_down
    !*************************************************************************************************









end module type_petsc_linear
