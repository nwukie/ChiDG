module type_chidg_vector_recv
#include <messenger.h>
    use mod_kinds,                   only: ik
    use type_ivector,                only: ivector_t
    use type_mesh,                   only: mesh_t
    use type_mesh_new,               only: mesh_new_t
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
        !type(mesh_t),               intent(inout)   :: mesh(:)
        type(mesh_new_t),           intent(inout)   :: mesh

        integer(ik)                 :: idom, iproc, icomm, ncomm, ndom_recv, ierr, loc
        integer(ik),    allocatable :: comm_procs_dom(:)
        type(ivector_t)             :: comm_procs
        logical                     :: not_in_list


        !
        ! Deallocate storage if necessary in case this is being called as a 
        ! reinitialization routine.
        !
        if (allocated(self%comm)) deallocate(self%comm)


        !
        ! Get processor ranks that we are receiving from: mesh
        !
        do idom = 1,mesh%ndomains()
            comm_procs_dom = mesh%domain(idom)%get_recv_procs()

            do iproc = 1,size(comm_procs_dom)
                ! See if proc is already in list
                loc = comm_procs%loc(comm_procs_dom(iproc))
                not_in_list = ( loc == 0 )

                ! If not, add to list
                if ( not_in_list ) call comm_procs%push_back(comm_procs_dom(iproc))
            end do ! iproc

        end do ! idom



        !
        ! Allocate recv communication for each processor sending data here.
        !
        ncomm = comm_procs%size()
        if (allocated(self%comm)) deallocate(self%comm)
        allocate(self%comm(ncomm), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Call initialization for each communicating process
        !
        do icomm = 1,comm_procs%size()
            iproc = comm_procs%at(icomm)
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
