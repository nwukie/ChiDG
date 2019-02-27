module operator_chidg_dot
#include <messenger.h>
#include "petsc/finclude/petscvec.h"
    use petscvec,               only: tVec, VecDot, VecGetLocalVector, VecCreate, VecSetType

    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO
    
    use type_chidg_vector,      only: chidg_vector_t
    use operator_domain_dot,    only: domain_dot
    use mod_chidg_mpi,          only: GROUP_MASTER
    use mpi_f08,                only: mpi_comm, MPI_AllReduce, MPI_Reduce, MPI_REAL8, MPI_SUM
    implicit none

    interface dot
        module procedure    dot_local, dot_comm
    end interface

contains



    !>  Compute vector-vector dot product from two chidg_vector_t types.
    !!  Computed from only the processor-local instances.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-----------------------------------------------------------------------------
    function dot_local(a,b) result(res)
        type(chidg_vector_t),    intent(in)  :: a
        type(chidg_vector_t),    intent(in)  :: b

        real(rk)    :: res
        integer(ik) :: idom, ielem

        res = ZERO

        if (a%petsc_vector_created) then
            call chidg_signal(FATAL,'dot_local: processor-local dot-product not yet implemented for petsc.')
        else

            ! Compute vector dot-product
            do idom = 1,size(a%dom)
                res = res + domain_dot(a%dom(idom),b%dom(idom))
            end do

        end if

    end function dot_local
    !******************************************************************************





    !>  Compute vector-vector dot product from two chidg_vector_t types. This
    !!  is computed across the processors within the communicator comm.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/23/2016
    !!
    !-----------------------------------------------------------------------------
    function dot_comm(a,b,comm) result(comm_dot)
        type(chidg_vector_t),   intent(in)  :: a
        type(chidg_vector_t),   intent(in)  :: b
        type(mpi_comm),         intent(in)  :: comm

        real(rk)    :: local_dot, comm_dot
        integer     :: ierr

        PetscErrorCode :: perr

        if (a%petsc_vector_created) then
            
            call VecDot(a%petsc_vector,b%petsc_vector,comm_dot,perr)
            if (perr /= 0) call chidg_signal(FATAL,'dot_comm: error calling petsc VecDot.')

        else

            ! Compute the local vector dot-product
            local_dot = dot_local(a,b)

            ! Reduce local dot-product values across processors, distribute result back to all
            call MPI_AllReduce(local_dot,comm_dot,1,MPI_REAL8,MPI_SUM,comm,ierr)

        end if

    end function dot_comm
    !******************************************************************************









end module operator_chidg_dot
