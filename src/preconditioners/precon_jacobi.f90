module precon_jacobi
#include <messenger.h>
#include "petsc/finclude/petscksp.h"
    use petscksp,               only: tPC, PCCreate, PCApply, PCDestroy, PCSetUp

    use mod_kinds,              only: rk, ik
    use mod_io,                 only: backend
    use mod_chidg_mpi,          only: ChiDG_COMM
    use type_preconditioner,    only: preconditioner_t
    use type_chidg_matrix,      only: chidg_matrix_t, chidg_matrix
    use type_chidg_vector,      only: chidg_vector_t
    use type_chidg_data,        only: chidg_data_t
    use mod_inv,                only: inv



    !> Block-Jacobi preconditioner
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!
    !-----------------------------------------------------------------------------
    type, extends(preconditioner_t) :: precon_jacobi_t

        type(chidg_matrix_t) :: D     ! inverse of block diagonal, (ndom,maxelems,ntime)

        PC      :: pc
        logical :: petsc_initialized = .false.

    contains
        procedure   :: init
        procedure   :: update
        procedure   :: apply
        procedure   :: restrict
        procedure   :: tear_down
    end type precon_jacobi_t
    !******************************************************************************




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
        class(precon_jacobi_t), intent(inout)   :: self
        type(chidg_data_t),     intent(in)      :: data

        PetscErrorCode :: perr
        
        select case (trim(backend))
            case('native')
                self%D = chidg_matrix(trim(backend))
                call self%D%init(data%mesh,mtype='Diagonal')

            case('petsc')

                call PCCreate(ChiDG_COMM%mpi_val,self%pc,perr)
                if (perr /= 0) call chidg_signal(FATAL,'precon_jacobi%init: error calling PCCreate.')
                call PCSetType(self%pc,PCPBJACOBI,perr)
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
    !!  @date   2/24/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update(self,A,b)
        class(precon_jacobi_t), intent(inout)   :: self
        type(chidg_matrix_t),    intent(in)      :: A
        type(chidg_vector_t),    intent(in)      :: b


        integer(ik)     :: idom, ielem, itime, diag
        PetscErrorCode  :: perr

        if (self%petsc_initialized) then
        !******  petsc  implementation  ******!
            call PCSetOperators(self%pc, A%wrapped_petsc_matrix%petsc_matrix, A%wrapped_petsc_matrix%petsc_matrix, perr)
            if (perr /= 0) call chidg_signal(FATAL,'precon_jacobi%update: error calling PCSetOperators.')
            call PCSetUp(self%pc, perr)
            if (perr /= 0) call chidg_signal(FATAL,'precon_jacobi%update: error calling PCSetUp.')


        else
        !******  ChiDG native implementation    ******!

            ! Copy the block diagonal from A to the preconditioner storage
            do idom = 1,size(A%dom)
                do ielem = 1,size(A%dom(idom)%lblks,1)
                    do itime = 1,size(A%dom(idom)%lblks,2)
                        diag = A%dom(idom)%lblks(ielem,itime)%get_diagonal()
                        self%D%dom(idom)%lblks(ielem,itime)%data_(1)%mat = A%dom(idom)%lblks(ielem,itime)%dmat(diag)
                    end do
                end do
            end do


            ! Replace the block diagonal D with Dinv
            do idom = 1,size(A%dom)
                do ielem = 1,size(A%dom(idom)%lblks,1)
                    do itime = 1,size(A%dom(idom)%lblks,2)
                        self%D%dom(idom)%lblks(ielem,itime)%data_(1)%mat = inv(self%D%dom(idom)%lblks(ielem,itime)%data_(1)%mat)
                    end do
                end do
            end do

        end if


        ! Update stamp
        self%stamp = A%stamp

    end subroutine update
    !***************************************************************************************








    !> Apply the preconditioner to the krylov vector 'v' and return preconditioned vector 'z'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    function apply(self,A,v,z_old) result(z)
        class(precon_jacobi_t), intent(inout)           :: self
        type(chidg_matrix_t),   intent(in)              :: A
        type(chidg_vector_t),   intent(in)              :: v
        type(chidg_vector_t),   intent(in), optional    :: z_old

        type(chidg_vector_t)    :: z
        integer(ik)             :: idom, ielem, itime, ndom
        real(rk),   allocatable :: mv(:)

        PetscErrorCode :: perr

        call self%timer%start()


        ! Allocate 'z'
        z = v
        call z%clear()


        if (self%petsc_initialized) then
        !******  petsc  implementation  ******!

            call PCApply(self%pc,v%wrapped_petsc_vector%petsc_vector,z%wrapped_petsc_vector%petsc_vector,perr)
            if (perr /= 0) call chidg_signal(FATAL,'precon_jacobi%apply: error calling PCApply.')
            z%petsc_needs_assembled = .true.


        else
        !******  chidg implementation  ******!

            ! Apply Block-Diagonal preconditioner
            do idom = 1,size(A%dom)
                do ielem = 1,size(A%dom(idom)%lblks,1)
                    do itime = 1,size(A%dom(idom)%lblks,2)
                        !z%dom(idom)%vecs(ielem)%vec = &
                        !    matmul(self%D%dom(idom)%lblks(ielem,itime)%data_(1)%mat,v%dom(idom)%vecs(ielem)%vec)
                        mv = matmul(self%D%dom(idom)%lblks(ielem,itime)%data_(1)%mat,v%dom(idom)%vecs(ielem)%gettime(itime))
                        call z%dom(idom)%vecs(ielem)%settime(itime, mv)
                    end do
                end do
            end do

        end if


        call self%timer%stop()

    end function apply
    !****************************************************************************************





    !>  Produce a restricted version of the current preconditioner.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/24/2017
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function restrict(self,nterms_r) result(restricted)
        class(precon_jacobi_t), intent(in)  :: self
        integer(ik),            intent(in)  :: nterms_r

        type(precon_jacobi_t) :: restricted

        restricted%D = self%D%restrict(nterms_r)
        restricted%initialized = .true.

    end function restrict
    !****************************************************************************************




    !>  Tear down. Deallocate, etc. 
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !-----------------------------------------------------------------------------
    subroutine tear_down(self)
        class(precon_jacobi_t), intent(inout)   :: self

        PetscErrorCode :: perr
        
        if (self%petsc_initialized) then
            call PCDestroy(self%pc, perr)
            if (perr /= 0) call chidg_signal(FATAL,'precon_jacobi%tear_down: error calling PCDestroy.')
            self%petsc_initialized = .false.

        else
            call self%D%release()

        end if

        ! Remove correspondense to matrix stamp
        self%stamp = 0

    end subroutine tear_down
    !***************************************************************************************










end module precon_jacobi
