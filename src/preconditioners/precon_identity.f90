module precon_identity
#include <messenger.h>
#include "petsc/finclude/petscksp.h"
    use petscksp,               only: tPC, PCCreate, PCApply, PCDestroy, PCSetUp

    use mod_kinds,              only: rk, ik
    use mod_io,                 only: backend
    use mod_constants,          only: DIAG, ONE
    use type_preconditioner,    only: preconditioner_t
    use type_chidg_data,        only: chidg_data_t
    use type_chidg_matrix,      only: chidg_matrix_t
    use type_chidg_vector


    !> Identity preconditioner
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !----------------------------------------------------------------------
    type, extends(preconditioner_t) :: precon_identity_t

        PC      :: pc
        logical :: petsc_initialized = .false.


    contains
        procedure   :: update
        procedure   :: apply

    end type precon_identity_t




contains


    !> Initialize preconditioner storage
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------
    subroutine init(self,data)
        class(precon_identity_t), intent(inout)   :: self
        type(chidg_data_t),     intent(in)      :: data

        PetscErrorCode :: perr
        
        select case (trim(backend))
            case('native')

            case('petsc')
                call PCCreate(ChiDG_COMM%mpi_val,self%pc,perr)
                if (perr /= 0) call chidg_signal(FATAL,'precon_jacobi%init: error calling PCCreate.')
                call PCSetType(self%pc,PCNONE,perr)
                if (perr /= 0) call chidg_signal(FATAL,'precon_jacobi%init: error calling PCSetType.')
                self%petsc_initialized = .true.

            case default
                call chidg_signal_one(FATAL,"precon_jacobi%init: invalid input for 'backend'.", trim(backend))

        end select

    end subroutine init
    !***************************************************************************************






    !> Compute the block diagonal inversion and store so it can be applied.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !---------------------------------------------------------------------
    subroutine update(self,A,b)
        class(precon_identity_t), intent(inout)   :: self
        type(chidg_matrix_t),    intent(in)      :: A
        type(chidg_vector_t),    intent(in)      :: b


        PetscErrorCode  :: perr

        if (self%petsc_initialized) then
        !******  petsc  implementation  ******!
            call PCSetOperators(self%pc, A%petsc_matrix, A%petsc_matrix, perr)
            if (perr /= 0) call chidg_signal(FATAL,'precon_jacobi%update: error calling PCSetOperators.')
            call PCSetUp(self%pc, perr)
            if (perr /= 0) call chidg_signal(FATAL,'precon_jacobi%update: error calling PCSetUp.')


        end if



        ! Update stamp
        call date_and_time(values=self%stamp)





        
    end subroutine update








    !> Apply the preconditioner to the krylov vector 'v' and return preconditioned vector 'z'
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !-------------------------------------------------------------------------
    function apply(self,A,v,z_old) result(z)
        class(precon_identity_t), intent(inout)   :: self
        type(chidg_matrix_t),     intent(in)      :: A
        type(chidg_vector_t),     intent(in)      :: v
        type(chidg_vector_t),     intent(in), optional :: z_old

        type(chidg_vector_t) :: z
        PetscErrorCode :: perr


        ! Identity preconditioner. Do nothing, copy incoming vector to outgoing vector without modification
        z = v


        if (self%petsc_initialized) then
        !******  petsc  implementation  ******!

            call PCApply(self%pc,v%petsc_vector,z%petsc_vector,perr)
            if (perr /= 0) call chidg_signal(FATAL,'precon_jacobi%apply: error calling PCApply.')
            z%from_operator         = .true.
            z%petsc_needs_assembled = .true.
        else

        end if


    end function apply
    !*************************************************************************













end module precon_identity
