module type_petsc_is_wrapper
#include <messenger.h>
#include "petsc/finclude/petscvec.h"
!#include "petscconf.h"
    use mod_kinds,  only: ik
    use petscis,    only: tIS, ISDestroy
    implicit none

    type, public :: petsc_is_wrapper_t
        IS :: petsc_is
    contains
        final :: destroy
    end type petsc_is_wrapper_t

contains


    !>  Deallocation of PETSc is. 
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/10/2019
    !!
    !------------------------------------------------------------------------------------
    subroutine destroy(self)
        type(petsc_is_wrapper_t),   intent(inout)   :: self

        integer(ik) :: ierr

        ! We assume, that whenever a petsc_is_wrapper_t is allocated, the petsc_is
        ! component is also created. So, we do not check that this is actually the case.
        call ISDestroy(self%petsc_is,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'petsc_is_wrapper%destroy: error calling ISDestroy.')


        ! NOTE: petsc sets vector objects to NULL after destroying them. The problem is that
        ! some petsc routines behave differently for NULL values.
        ! https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecisCreateToAll.html
        !
        ! So, we reinitialize the petsc internal
        ! representation here to PETSC_FORTRAN_TYPE_INITIALIZE, which is their default value
        !   See: src/vec/f90-mod/petscvec.h
        !
        ! Also, the PETSC_FORTRAN_TYPE_INITIALIZE macro is defined as '= -2' so there shouldn't
        ! be any '=' below or the compiler will error.
        !
        self%petsc_is%v      PETSC_FORTRAN_TYPE_INITIALIZE


    end subroutine destroy
    !*************************************************************************************


end module type_petsc_is_wrapper
