module type_chidg_vector_send
#include <messenger.h>
    use mod_kinds,                      only: ik
    use type_ivector,                   only: ivector_t
    use type_mesh_new,                  only: mesh_new_t
    use type_chidg_vector_send_comm,    only: chidg_vector_send_comm_t
    use mpi_f08,                        only: MPI_Request, MPI_STATUSES_IGNORE, MPI_Waitall
    implicit none




    !>  Container for storing information about what parts of chidg_vector to send to other
    !!  processors.
    !!
    !!  For each processor that we are sending data to, a chidg_vector_send_comm_t instance exists:
    !!      self%comm(icomm)
    !!
    !!  The chidg_vector_send_comm_t instance contains all the information about what
    !!  data we are sending to a specific processor.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    type, public :: chidg_vector_send_t

        type(chidg_vector_send_comm_t),  allocatable :: comm(:)
        type(MPI_Request),               allocatable :: isend_handles(:)

    contains

        procedure, public   :: init
        procedure, public   :: init_wait

    end type chidg_vector_send_t
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
        class(chidg_vector_send_t),  intent(inout)   :: self
        type(mesh_new_t),            intent(in)      :: mesh

        integer(ik) :: idom, iproc, nprocs_send, ierr, loc, nsends
        integer(ik), allocatable    :: comm_procs_dom(:)
        type(ivector_t)             :: comm_procs
        logical                     :: already_added


        !
        ! Deallocate storage if necessary in case this is being called as a 
        ! reinitialization routine.
        !
        if (allocated(self%comm)) deallocate(self%comm)
        if (allocated(self%isend_handles)) deallocate(self%isend_handles)


        !
        ! Detect processors that we need to send to, from all mesh instances.
        !
        do idom = 1,mesh%ndomains()

            ! Get comm_procs from current domain
            comm_procs_dom = mesh%domain(idom)%get_send_procs() 
            
            ! Loop through and add them to list if not already added
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
        ! Accumulate total number of sends and allocate storage for mpi_requests for each send. 
        ! That way, chidg_vector%comm_wait can wait on them to complete.
        !
        nsends = 0
        do iproc = 1,comm_procs%size()
            nsends = nsends + self%comm(iproc)%nsends()
        end do ! iproc

        allocate(self%isend_handles(nsends), stat=ierr)
        if (ierr /= 0) call AllocationError

    end subroutine init
    !****************************************************************************************







    !>  Wait on any outstanding mpi_requests that were sent during 
    !!  initiailization, self%init.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/13/2017
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine init_wait(self)
        class(chidg_vector_send_t),  intent(inout)   :: self
        
        integer(ik) :: icomm, nrequests, ierr

        
        do icomm = 1,size(self%comm)

            nrequests = self%comm(icomm)%initialization_requests%size()
            if (nrequests > 0) then
                call MPI_Waitall(nrequests, self%comm(icomm)%initialization_requests%data, MPI_STATUSES_IGNORE, ierr)
            end if

        end do

    end subroutine init_wait
    !***************************************************************************************





end module type_chidg_vector_send
