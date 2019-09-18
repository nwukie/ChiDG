module type_petsc_matrix_wrapper
#include <messenger.h>
#include "petsc/finclude/petscmat.h"
    use mod_kinds,  only: ik
    use petscmat,   only: tMat, MatDestroy
    implicit none

    type, public :: petsc_matrix_wrapper_t
        Mat     :: petsc_matrix
    contains
        final :: destroy
    end type petsc_matrix_wrapper_t

contains


    !>  Deallocation of PETSc Vector. 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2019
    !!
    !------------------------------------------------------------------------------------
    subroutine destroy(self)
        type(petsc_matrix_wrapper_t),   intent(inout)   :: self

        integer(ik) :: ierr

        ! We assume, that whenever a petsc_vector_wrapper_t is allocated, the petsc_vector 
        ! component is also created. So, we do not check that this is actually the case.
        call MatDestroy(self%petsc_matrix,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'petsc_matrix_wrapper%destroy: error calling MatDestroy.')


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
        self%petsc_matrix%v      PETSC_FORTRAN_TYPE_INITIALIZE


    end subroutine destroy
    !*************************************************************************************


end module type_petsc_matrix_wrapper
