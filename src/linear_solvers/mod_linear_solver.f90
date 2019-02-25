module mod_linear_solver
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use type_linear_solver, only: linear_solver_t
    use type_dict,          only: dict_t

    use type_fgmres,                only: fgmres_t
    use type_petsc_linear,          only: petsc_linear_t
    use type_fgmres_cgs_mg_correct, only: fgmres_cgs_mg_correct_t
    implicit none

    type(fgmres_t)                  :: FGMRES
    type(petsc_linear_t)            :: PETSC_LINEAR
    type(fgmres_cgs_mg_correct_t)   :: FGMRES_CGS_MG_CORRECT

contains

    !>  Factory method for creating matrixsolver objects
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/17/2016
    !!
    !!  @param[in]      mstring     Character string used to select the appropriate matrixsolver for allocation
    !!  @param[inout]   msolver     matrixsolver_t that will be allocated to a concrete type.
    !!
    !---------------------------------------------------------------------------------------------------------------------
    subroutine create_linear_solver(lstring,lsolver,options)
        character(*),                        intent(in)      :: lstring
        class(linear_solver_t), allocatable, intent(inout)   :: lsolver
        type(dict_t), optional,              intent(inout)   :: options

        integer(ik) :: ierr


        select case (trim(lstring))
            case ('fgmres','FGMRES')
                allocate(lsolver, source=FGMRES, stat=ierr)

            case ('petsc','PETSC')
                allocate(lsolver, source=PETSC_LINEAR, stat=ierr)

            case ('fgmres_cgs_mg_correct', 'FGMRES_CGS_MG_CORRECT')
                allocate(lsolver, source=FGMRES_CGS_MG_CORRECT, stat=ierr)

            case default
                call chidg_signal(FATAL,"create_matrixsolver: matrix solver string did not match any valid type")

        end select
        if (ierr /= 0) call AllocationError


        ! Call options initialization if present
        if (present(options)) call lsolver%set(options)

        
        ! Make sure the solver was allocated
        if (.not. allocated(lsolver)) call chidg_signal(FATAL,"create_matrixsolver: solver was not allocated. Check that the desired solver was registered and instantiated in the mod_matrixsolver module")

    end subroutine create_linear_solver
    !*********************************************************************************************************************




end module mod_linear_solver
