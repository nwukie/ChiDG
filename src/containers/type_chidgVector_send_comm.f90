module type_chidgVector_send_comm
#include <messenger.h>
    use mod_kinds,          only: ik
    use type_mesh,          only: mesh_t
    use type_ivector,       only: ivector_t
    use mod_chidg_mpi,      only: ChiDG_COMM
    use mpi_f08,            only: MPI_Request, MPI_INTEGER4, MPI_ISend
    implicit none






    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/1/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    type, public :: chidgVector_send_comm_t

        integer(ik)                     :: proc
        type(ivector_t)                 :: dom_send         !< Vector of domain indices in the mesh that have elems to be sent.
        type(ivector_t),    allocatable :: elems_send(:)    !< For each domain with info to be sent, a 
                                                            !< vector that contains the indices of elements to be sent.
    contains

        procedure,  public  :: init     !< Initialize the send info going to a particular process
        procedure,  public  :: nsends   !< Return the number of sends going to the process

    end type chidgVector_send_comm_t
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
        class(chidgVector_send_comm_t), intent(inout)   :: self
        type(mesh_t),                   intent(in)      :: mesh(:)
        integer(ik),                    intent(in)      :: proc

        integer(ik)                 :: idom, ielem, iface, idom_send, ndom_send, ierr, loc, neighbor_proc, ielem_send
        integer(ik),    allocatable :: comm_procs_dom(:)
        logical                     :: already_added, proc_has_domain, proc_is_neighbor
        type(mpi_request)           :: null_request


        !
        ! Set processor being sent to
        !
        self%proc = proc


        !
        ! For the current proc that we are sending stuff to, detect how many domains are being sent
        !
        do idom = 1,size(mesh)

            ! Get comm procs for domain
            comm_procs_dom = mesh(idom)%get_comm_procs()

            ! Does proc we are sending to have comm with the current domain
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
            call self%dom_send%clear()
            deallocate(self%elems_send)
        end if
        allocate(self%elems_send(ndom_send), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Initialize element send indices for each domain to be sent
        !
        do idom_send = 1,self%dom_send%size()

            idom = self%dom_send%at(idom_send)

            do ielem = 1,mesh(idom)%nelem
                do iface = 1,size(mesh(idom)%faces,2)


                    ! Check if element has comm with neighbor at iface

                    neighbor_proc = mesh(idom)%faces(ielem,iface)%ineighbor_proc
                    proc_is_neighbor = ( proc == neighbor_proc )


                    ! If element should be sent, add to list
                    if (proc_is_neighbor) then
                        ! Check if element was already added to send
                        loc = self%elems_send(idom_send)%loc(ielem)
                        already_added = (loc /= 0)

                        ! Add to send list if not already there
                        if ( .not. already_added ) then
                            ! Add to send list
                            call self%elems_send(idom_send)%push_back(ielem)

                        end if
                    end if


                end do ! iface
            end do ! ielem
        end do ! idom









        !
        ! Initialize element send indices for each domain to be sent
        !
        do idom_send = 1,self%dom_send%size()

            idom = self%dom_send%at(idom_send)

            ! Communicate domain index
            call MPI_ISend(mesh(idom)%idomain_g, 1, MPI_INTEGER4, self%proc, 0, ChiDG_COMM, null_request, ierr)

            ! Communicate number of elements being sent
            call MPI_ISend(self%elems_send(idom_send)%size_, 1, MPI_INTEGER4, self%proc, 0, ChiDG_COMM, null_request, ierr)




            do ielem_send = 1,self%elems_send(idom_send)%size()
                ielem = self%elems_send(idom_send)%at(ielem_send)

                ! Communicate element index
                call MPI_ISend(mesh(idom)%elems(ielem)%ielement_g, 1, MPI_INTEGER4, self%proc, 0, ChiDG_COMM, null_request, ierr)
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
        class(chidgVector_send_comm_t), intent(in)  :: self

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





end module type_chidgVector_send_comm
