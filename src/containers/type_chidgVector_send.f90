module type_chidgVector_send
#include <messenger.h>
    use mod_kinds,                  only: ik
    use type_ivector,               only: ivector_t
    use type_mesh,                  only: mesh_t
    use type_chidgVector_send_comm, only: chidgVector_send_comm_t
    use mpi_f08,                    only: MPI_Request
    implicit none


    !>  Container for storing information about what parts of chidgVector to send to other
    !!  processors.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    type, public :: chidgVector_send_t

        type(chidgVector_send_comm_t),  allocatable :: comm(:)
        type(MPI_Request),              allocatable :: isend_handles(:)

    contains

        procedure, public :: init

    end type chidgVector_send_t
    !*****************************************************************************************






contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine init(self,mesh)
        class(chidgVector_send_t),  intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh(:)

        integer(ik) :: idom, iproc, nprocs_send, ierr, loc, nsends
        integer(ik), allocatable    :: comm_procs_dom(:)
        type(ivector_t)             :: comm_procs
        logical                     :: already_added


        !
        ! Detect number of processors that we need to send to
        !
        do idom = 1,size(mesh)

            comm_procs_dom = mesh(idom)%get_comm_procs() 
            
            do iproc = 1,size(comm_procs_dom)
                ! check if proc was already added to list from another domain
                loc = comm_procs%loc(comm_procs_dom(iproc))
                already_added = (loc /= 0)
                if ( .not. already_added ) call comm_procs%push_back(comm_procs_dom(iproc))
            end do ! iproc

        end do ! idom


        !
        ! Get total number of procs we are sending to and allocate send info for each
        !
        nprocs_send = comm_procs%size()
        allocate(self%comm(nprocs_send), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Call initialization for each proc that we are sending to
        !
        do iproc = 1,comm_procs%size()
            call self%comm(iproc)%init(mesh,comm_procs%at(iproc))
        end do !iproc



        !
        ! Accumulate total number of sends and allocate storage for mpi_requests for each send
        !
        nsends = 0
        do iproc = 1,comm_procs%size()
            nsends = nsends + self%comm(iproc)%nsends()
        end do ! iproc

        allocate(self%isend_handles(nsends), stat=ierr)
        if (ierr /= 0) call AllocationError

    end subroutine init
    !****************************************************************************************









end module type_chidgVector_send
