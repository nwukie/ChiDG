module type_chidg_vector_send_comm
#include <messenger.h>
    use mod_kinds,                  only: ik
    use mod_constants,              only: INTERIOR, CHIMERA
    use type_mesh,              only: mesh_t
    use type_ivector,               only: ivector_t
    use mod_chidg_mpi,              only: ChiDG_COMM
    use type_mpi_request_vector,    only: mpi_request_vector_t
    use mpi_f08,                    only: MPI_Request, MPI_INTEGER4, MPI_ISend
    implicit none






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    type, public :: chidg_vector_send_comm_t

        integer(ik)                     :: proc
        type(ivector_t)                 :: dom_send         ! Vector of domain indices in the mesh that have elems to be sent.
        type(ivector_t),    allocatable :: elems_send(:)    ! For each domain with info to be sent, a 
                                                            ! vector that contains the indices of elements to be sent.

        type(mpi_request_vector_t)      :: initialization_requests

    contains

        procedure,  public  :: init     ! Initialize the send info going to a particular process
        procedure,  public  :: nsends   ! Return the number of sends going to the process

    end type chidg_vector_send_comm_t
    !************************************************************************************












contains





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------
    subroutine init(self,mesh,proc)
        class(chidg_vector_send_comm_t),    intent(inout)   :: self
        type(mesh_t),                       intent(in)      :: mesh
        integer(ik),                        intent(in)      :: proc

        integer(ik)                 :: idom, ielem, iface, idom_send, ndom_send, ierr, &
                                       loc, neighbor_proc, ielem_send, ChiID, idonor, receiver_proc
        integer(ik),    allocatable :: comm_procs_dom(:)
        logical                     :: already_added, proc_has_domain, send_element, &
                                       has_neighbor, is_chimera
        type(mpi_request)           :: request, request1, request2, request3, request4, request5


        !
        ! Set processor being sent to
        !
        self%proc = proc


        !
        ! For the current proc that we are sending stuff to, detect how many domains are being sent
        !
        do idom = 1,mesh%ndomains()

            ! Get send procs for domain
            comm_procs_dom = mesh%domain(idom)%get_send_procs()

            ! Does current domain have send comm with proc we are sending to here
            proc_has_domain = any( proc == comm_procs_dom )

            ! If so, add domain to list of domains being sent to proc
            if (proc_has_domain) then
                call self%dom_send%push_back(idom)
            end if

        end do ! idom



        !
        ! Allocate number of domains to send to proc
        !
        ndom_send = self%dom_send%size()
        if (allocated(self%elems_send)) then
            deallocate(self%elems_send)
        end if
        allocate(self%elems_send(ndom_send), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Initialize element send indices for each domain to be sent
        !
        do idom_send = 1,self%dom_send%size()


            idom = self%dom_send%at(idom_send)


            !
            ! Register elements to send to off-processor neighbors
            !
            do ielem = 1,mesh%domain(idom)%nelem
                do iface = 1,size(mesh%domain(idom)%faces,2)

                    has_neighbor = ( mesh%domain(idom)%faces(ielem,iface)%ftype == INTERIOR )
                    if ( has_neighbor ) then
                        !
                        ! Get neighbor processor rank
                        !
                        neighbor_proc = mesh%domain(idom)%faces(ielem,iface)%ineighbor_proc
                        send_element = ( proc == neighbor_proc )

                        !
                        ! If element should be sent, add to list
                        !
                        if (send_element) then
                            ! Check if element was already added to send
                            loc = self%elems_send(idom_send)%loc(ielem)
                            already_added = (loc /= 0)

                            ! Add to send list if not already there
                            if ( .not. already_added ) call self%elems_send(idom_send)%push_back(ielem)
                        end if

                    end if

                end do ! iface
            end do ! ielem




            !
            ! Register elements to send to off-processor chimera receivers 
            !
            do idonor = 1,mesh%domain(idom)%chimera%send%ndonors()

                ! Get proc of receiver
                receiver_proc = mesh%domain(idom)%chimera%send%receiver_proc%at(idonor)
                send_element  = (proc == receiver_proc)

                !
                ! If element should be sent, add to list
                !
                if (send_element) then
                    ielem = mesh%domain(idom)%chimera%send%donor_element_l%at(idonor)
                    ! Check if element was already added to send
                    loc = self%elems_send(idom_send)%loc(ielem)
                    already_added = (loc /= 0)

                    ! Add to send list if not already there
                    if ( .not. already_added ) call self%elems_send(idom_send)%push_back(ielem)
                end if


            end do !idonor


        end do ! idom









        !
        ! Initialize element send indices for each domain to be sent
        !

        ! Communicate number of domains being sent to proc. This is recv'd by chidg_vector_recv_comm
        call MPI_ISend(self%dom_send%size_, 1, MPI_INTEGER4, self%proc, 0, ChiDG_COMM, request, ierr)
        call self%initialization_requests%push_back(request)


        ! These send's are recv'd by blockvector%init_recv
        do idom_send = 1,self%dom_send%size()

            idom = self%dom_send%at(idom_send)

            ! Communicate domain indices
            call MPI_ISend(mesh%domain(idom)%idomain_g, 1, MPI_INTEGER4, self%proc, 0, ChiDG_COMM, request1, ierr)
            call MPI_ISend(mesh%domain(idom)%idomain_l, 1, MPI_INTEGER4, self%proc, 0, ChiDG_COMM, request2, ierr)

            call self%initialization_requests%push_back(request1)
            call self%initialization_requests%push_back(request2)

            ! Communicate number of elements from the domain being sent
            call MPI_ISend(self%elems_send(idom_send)%size_, 1, MPI_INTEGER4, self%proc, 0, ChiDG_COMM, request, ierr)
            call self%initialization_requests%push_back(request)

            ! Send each element
            do ielem_send = 1,self%elems_send(idom_send)%size()
                ielem = self%elems_send(idom_send)%at(ielem_send)

                call MPI_ISend(mesh%domain(idom)%elems(ielem)%ielement_g, 1, MPI_INTEGER4, self%proc, 0, ChiDG_COMM, request1, ierr)
                call MPI_ISend(mesh%domain(idom)%elems(ielem)%ielement_l, 1, MPI_INTEGER4, self%proc, 0, ChiDG_COMM, request2, ierr)
                call MPI_ISend(mesh%domain(idom)%elems(ielem)%nterms_s,   1, MPI_INTEGER4, self%proc, 0, ChiDG_COMM, request3, ierr)
                call MPI_ISend(mesh%domain(idom)%elems(ielem)%neqns,      1, MPI_INTEGER4, self%proc, 0, ChiDG_COMM, request4, ierr)
                call MPI_ISend(mesh%domain(idom)%elems(ielem)%ntime,      1, MPI_INTEGER4, self%proc, 0, ChiDG_COMM, request5, ierr)

                call self%initialization_requests%push_back(request1)
                call self%initialization_requests%push_back(request2)
                call self%initialization_requests%push_back(request3)
                call self%initialization_requests%push_back(request4)
                call self%initialization_requests%push_back(request5)

            end do ! ielem_send

        end do ! idom_send




    end subroutine init
    !*************************************************************************************











    !>  Return the number of sends going to the process associated with the current container
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    function nsends(self) result(nsend)
        class(chidg_vector_send_comm_t), intent(in)  :: self

        integer(ik) :: nsend, idom_send, ielem_send

        !
        ! Loop through domains
        !
        nsend = 0
        do idom_send = 1,self%dom_send%size()

            ! Accumulate number of elements being sent from each sending domain
            nsend = nsend + self%elems_send(idom_send)%size()

        end do

    end function nsends
    !****************************************************************************************





end module type_chidg_vector_send_comm
