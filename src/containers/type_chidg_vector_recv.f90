module type_chidg_vector_recv
#include <messenger.h>
    use mod_kinds,                   only: ik
    use type_ivector,                only: ivector_t
    use type_mesh,                   only: mesh_t
    use type_chidg_vector_recv_comm, only: chidg_vector_recv_comm_t
    implicit none



    !>  Container for receiving vector data from other chidg_vectors. Holds a container for
    !!  each process that is sending data here.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    type, public :: chidg_vector_recv_t

        type(chidg_vector_recv_comm_t),   allocatable :: comm(:)


    contains

        procedure, public :: init
        procedure, public :: clear

    end type chidg_vector_recv_t
    !*****************************************************************************************





contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self,mesh)
        class(chidg_vector_recv_t), intent(inout)   :: self
        type(mesh_t),               intent(inout)   :: mesh

        integer(ik)                 :: idom, iproc, icomm, ncomm, ierr
        integer(ik),    allocatable :: comm_procs_array(:)


        !
        ! Deallocate storage if necessary in case this is being called as a 
        ! reinitialization routine.
        !
        if (allocated(self%comm)) deallocate(self%comm)


        !
        ! Get processor ranks that we are receiving from: mesh
        !
        comm_procs_array = mesh%get_recv_procs()


        !
        ! Allocate recv communication for each processor sending data here.
        !
        ncomm = size(comm_procs_array)
        if (allocated(self%comm)) deallocate(self%comm)
        allocate(self%comm(ncomm), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Call initialization for each communicating process
        !
        do icomm = 1,ncomm
            iproc = comm_procs_array(icomm)
            call self%comm(icomm)%init(mesh,iproc,icomm)
        end do



    end subroutine init
    !*******************************************************************************************







    !>  Clear the received data, but do NOT deallocated.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/7/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(chidg_vector_recv_t),  intent(inout)   :: self
        
        integer(ik) :: icomm

    
        do icomm = 1,size(self%comm)
            call self%comm(icomm)%clear()
        end do !icomm


    end subroutine clear
    !********************************************************************************************



end module type_chidg_vector_recv
