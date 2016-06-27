module operator_chidg_dot
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: ZERO
    
    use type_chidgVector,   only: chidgVector_t
    use operator_block_dot, only: block_dot
    use mod_chidg_mpi,      only: GROUP_MASTER
    use mpi_f08,            only: mpi_comm, MPI_Reduce, MPI_REAL8, MPI_SUM
    implicit none

    interface dot
        module procedure    dot_local, dot_comm
    end interface

contains



    !>  Compute vector-vector dot product from two chidgVector_t types.
    !!  Computed from only the processor-local instances.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-----------------------------------------------------------------------------
    function dot_local(a,b) result(res)
        type(chidgVector_t),    intent(in)  :: a
        type(chidgVector_t),    intent(in)  :: b

        real(rk)    :: res
        integer(ik) :: idom, ielem

        res = ZERO

        ! Compute vector dot-product
        do idom = 1,size(a%dom)
            
            res = res + block_dot(a%dom(idom),b%dom(idom))

        end do

    end function dot_local
    !******************************************************************************





    !>  Compute vector-vector dot product from two chidgVector_t types. This
    !!  is computed across the processors within the communicator comm.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/23/2016
    !!
    !-----------------------------------------------------------------------------
    function dot_comm(a,b,comm) result(comm_dot)
        type(chidgVector_t),    intent(in)  :: a
        type(chidgVector_t),    intent(in)  :: b
        type(mpi_comm),         intent(in)  :: comm

        real(rk)    :: local_dot, comm_dot
        integer(ik) :: ierr


        ! Compute the local vector dot-product
        local_dot = dot_local(a,b)


        ! Reduce local dot-product values across processors
        call MPI_Reduce(local_dot,comm_dot,1,MPI_REAL8,MPI_SUM,GROUP_MASTER,comm,ierr)


        ! Broadcast result to the group
        call MPI_BCast(comm_dot,1,MPI_REAL8,GROUP_MASTER,comm,ierr)


    end function dot_comm
    !******************************************************************************









end module operator_chidg_dot
