module precon_petsc
#include <messenger.h>
#include "petsc/finclude/petscksp.h"
    use petscksp,               only: tPC, PCCreate, PCApply, PCDestroy, PCSetUp

    use mod_kinds,              only: rk, ik
    use mod_constants,          only: DIAG, XI_MIN, ETA_MIN, ZETA_MIN, XI_MAX, ETA_MAX, ZETA_MAX, ONE
    use mod_inv,                only: inv
    use mod_chidg_mpi,          only: IRANK
    use mod_io,                 only: verbosity, backend

    use type_preconditioner,    only: preconditioner_t
    use type_chidg_data,        only: chidg_data_t
    use type_chidg_matrix,      only: chidg_matrix_t, chidg_matrix
    use type_chidg_vector,      only: chidg_vector_t
    implicit none





    !>  Use petsc preconditioner that can be configured from the command line.
    !!
    !!  Command line examples:
    !!      -pc_type asm -sub_pc_type ilu -sub_pc_factor_levels 1 -sub_ksp_type preonly
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/27/2019
    !!
    !!
    !-------------------------------------------------------------------------------------------
    type, extends(preconditioner_t) :: precon_petsc_t

        PC      :: pc
        logical :: petsc_initialized = .false.

    contains
    
        procedure   :: init
        procedure   :: update
        procedure   :: apply
        procedure   :: tear_down
        procedure   :: restrict

    end type precon_petsc_t
    !*******************************************************************************************




contains



    !>  Initialize the ILU0 preconditioner. This is for allocating storage. In this case, 
    !!  we allocate a Lower-Diagonal block matrix for storing the LU decomposition.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!  @param[inout]   domain  domain_t instance containing a mesh component used to 
    !!                          initialize the block matrix
    !!
    !-------------------------------------------------------------------------------------------
    subroutine init(self,data)
        class(precon_petsc_t),  intent(inout)   :: self
        type(chidg_data_t),     intent(in)      :: data

        PetscErrorCode :: perr
        
        select case (trim(backend))
            case('native')

                call chidg_signal(FATAL,"precon_petsc%init: tried to use petsc preconditioner with native containers. Must use backend = 'petsc'.")

            case('petsc')

                call PCCreate(ChiDG_COMM%mpi_val,self%pc,perr)
                if (perr /= 0) call chidg_signal(FATAL,'precon_petsc%init: error calling PCCreate.')

                call PCSetOperators(self%pc,data%sdata%lhs%wrapped_petsc_matrix%petsc_matrix,data%sdata%lhs%wrapped_petsc_matrix%petsc_matrix,perr)
                if (perr /= 0) call chidg_signal(FATAL,'petsc_petsc%init: error calling PCSetOperators.')

                call PCSetFromOptions(self%pc,perr)
                if (perr /= 0) call chidg_signal(FATAL,'petsc_petsc%init: error calling PCSetFromOptions.')

                self%petsc_initialized = .true.

            case default
                call chidg_signal_one(FATAL,"precon_petsc%init: invalid input for 'backend'.", trim(backend))

        end select

        self%initialized = .true.

    end subroutine init
    !*******************************************************************************************








    !>  Compute the block diagonal inversion and store so it can be applied.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!  Updated for spectral time integrators
    !!
    !!  @author Mayank Sharma + Nathan Wukie
    !!  @date   10/09/2017
    !!
    !-------------------------------------------------------------------------------------------
    subroutine update(self,A,b)
        class(precon_petsc_t),  intent(inout)   :: self
        type(chidg_matrix_t),   intent(in)      :: A
        type(chidg_vector_t),   intent(in)      :: b


        integer(ik) :: idom, ielem, itime, idiagA, idiagLD, irow, icol, &
                       eparent_l, ilowerA, ilowerLD, itranspose, dparent_g_lower, &
                       eparent_g_lower

        PetscErrorCode  :: perr

        call write_line(' Updating preconditioner...', io_proc=GLOBAL_MASTER, silence=(verbosity<5))

        call PCSetOperators(self%pc, A%wrapped_petsc_matrix%petsc_matrix, A%wrapped_petsc_matrix%petsc_matrix, perr)
        if (perr /= 0) call chidg_signal(FATAL,'precon_petsc%update: error calling PCSetOperators.')
        call PCSetUp(self%pc, perr)
        if (perr /= 0) call chidg_signal(FATAL,'precon_petsc%update: error calling PCSetUp.')


        call write_line(' Done updating preconditioner...', io_proc=GLOBAL_MASTER, silence=(verbosity<5))

        ! Update stamp
        self%stamp = A%stamp


    end subroutine update
    !*******************************************************************************************








    !> Apply the preconditioner to the krylov vector 'v' and return preconditioned vector 'z'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!  Updated for spectral time integrators
    !!
    !!  @author Mayank Sharma + Nathan Wukie
    !!  @date   10/09/2017
    !!
    !-------------------------------------------------------------------------------------------
    function apply(self,A,v,z_old) result(z)
        class(precon_petsc_t),   intent(inout)   :: self
        type(chidg_matrix_t),    intent(in)      :: A
        type(chidg_vector_t),    intent(in)      :: v
        type(chidg_vector_t),    intent(in), optional :: z_old

        type(chidg_vector_t)         :: z

        PetscErrorCode          :: perr

        call self%timer%start()

        ! Initialize z for preconditioning
        z = v

        call PCApply(self%pc,v%wrapped_petsc_vector%petsc_vector,z%wrapped_petsc_vector%petsc_vector,perr)
        if (perr /= 0) call chidg_signal(FATAL,'precon_petsc%apply: error calling PCApply.')

        call self%timer%stop()

    end function apply
    !-----------------------------------------------------------------------------------------








    !>  Produce a restricted version of the current preconditioner.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/24/2017
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function restrict(self,nterms_r) result(restricted)
        class(precon_petsc_t),   intent(in)  :: self
        integer(ik),            intent(in)  :: nterms_r

        type(precon_petsc_t)    :: restricted

        call chidg_signal(FATAL,'precon_petsc%restrict: not yet implemented for petsc preconditioners.')

    end function restrict
    !****************************************************************************************





    !>  Tear down. Deallocate, etc. 
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !-----------------------------------------------------------------------------
    subroutine tear_down(self)
        class(precon_petsc_t), intent(inout)   :: self

        PetscErrorCode :: perr
        
        if (self%petsc_initialized) then
            call PCDestroy(self%pc, perr)
            if (perr /= 0) call chidg_signal(FATAL,'precon_petsc%tear_down: error calling PCDestroy.')
            self%petsc_initialized = .false.
        end if

        ! Remove correspondense to matrix stamp
        self%stamp = 0

    end subroutine tear_down
    !***************************************************************************************





end module precon_petsc
