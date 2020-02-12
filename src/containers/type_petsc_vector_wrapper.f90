module type_petsc_vector_wrapper
#include <messenger.h>
#include "petsc/finclude/petscvec.h"
    use mod_kinds,  only: ik
    use petscvec,   only: tVec, VecDestroy
    implicit none

    type, public :: petsc_vector_wrapper_t
        Vec :: petsc_vector
    contains
        final :: destroy
    end type petsc_vector_wrapper_t

contains


    !>  Deallocation of PETSc Vector. 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2019
    !!
    !------------------------------------------------------------------------------------
    subroutine destroy(self)
        type(petsc_vector_wrapper_t),   intent(inout)   :: self

        PetscErrorCode :: ierr

        ! We assume, that whenever a petsc_vector_wrapper_t is allocated, the petsc_vector 
        ! component is also created. So, we do not check that this is actually the case.
        call VecDestroy(self%petsc_vector,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'petsc_vector_wrapper%destroy: error calling VecDestroy.')

        ! NOTE: petsc sets vector objects to NULL after destroying them. The problem is that
        ! VecScatterCreateToAll won't create a new vector for petsc_vector_recv if it is NULL
        ! because it thinks you passed NULL in as a flag. 
        ! https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecScatterCreateToAll.html
        !
        ! So, we reinitialize the petsc internal
        ! representation here to PETSC_FORTRAN_TYPE_INITIALIZE, which is their default value
        !   See: src/vec/f90-mod/petscvec.h
        !
        ! Also, the PETSC_FORTRAN_TYPE_INITIALIZE macro is defined as '= -2' so there shouldn't
        ! be any '=' below or the compiler will error.
        !
        self%petsc_vector%v      PETSC_FORTRAN_TYPE_INITIALIZE

    end subroutine destroy
    !*************************************************************************************


end module type_petsc_vector_wrapper
