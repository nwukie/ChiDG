module type_chidgVector_recv_comm
#include <messenger.h>
    use mod_kinds,          only: ik
    use type_mesh,          only: mesh_t
    use type_ivector,       only: ivector_t
    use type_blockvector,   only: blockvector_t

    !test
    use mod_chidg_mpi,      only: ChiDG_COMM
    use mpi_f08,            only: MPI_Recv, MPI_INTEGER4, MPI_STATUS_IGNORE
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !!
    !----------------------------------------------------------------------------
    type, public :: chidgVector_recv_comm_t

        integer(ik)                         :: proc
        type(blockvector_t),    allocatable :: dom(:)


        ! Mappings for where to store incoming data
        type(ivector_t)                 :: dom_store
        type(ivector_t),    allocatable :: elem_store(:)


    contains

        procedure,  public  :: init
        procedure,  public  :: clear

    end type chidgVector_recv_comm_t
    !****************************************************************************








contains




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine init(self,mesh,iproc,icomm)
        class(chidgVector_recv_comm_t), intent(inout)   :: self
        type(mesh_t),                   intent(inout)   :: mesh(:)
        integer(ik),                    intent(in)      :: iproc
        integer(ik),                    intent(in)      :: icomm

        type(ivector_t)             :: dom_recv
        integer(ik)                 :: idom, idom_recv, ndom_recv, ierr
        integer(ik)                 :: idomain_g, ielement_g, dom_store, idom_loop, ielem_loop, ielem_recv, nelem_recv
        integer(ik),    allocatable :: comm_procs_dom(:)
        logical                     :: proc_has_domain, domain_found, element_found

        !
        ! Set processor being received from
        !
        self%proc = iproc


        !
        ! Compute number of domains being received from proc
        !
        do idom = 1,size(mesh)

            ! Get comm procs for domain
            comm_procs_dom = mesh(idom)%get_comm_procs()

            ! Is proc communicating any of the current domain
            proc_has_domain = any( iproc == comm_procs_dom )

            ! If so, add domain to list of domains being received from current_proc
            if (proc_has_domain) then
                call dom_recv%push_back(idom)
            end if

        end do ! idom



        !
        ! Allocate number of domains to recv from proc
        !
        ndom_recv = dom_recv%size()
        if (allocated(self%dom)) deallocate(self%dom)
        allocate(self%dom(ndom_recv), self%elem_store(ndom_recv), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Call initialization for each domain to be received
        !
        do idom_recv = 1,dom_recv%size()
            call self%dom(idom_recv)%init(mesh(dom_recv%at(idom_recv)),iproc,icomm,idom_recv)
        end do ! idom




        !
        ! Get information from send processor about ordering of sends
        !
        ! Because: on this processor, we might detect the off-processor neighbors to recv in a 
        !          different order than they were detected on the neighbor processor to be sent.
        !          So, their storage order could be different on the processors. So, we communicate a bit,
        !          and store some integer maps to get the correct domain and element indices to store 
        !          data being recv'd.
        !
        do idom_recv = 1,dom_recv%size()

            call MPI_Recv(idomain_g,  1, MPI_INTEGER4, self%proc, 0, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
            call MPI_Recv(nelem_recv, 1, MPI_INTEGER4, self%proc, 0, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)


            ! Find domain store index
            domain_found = .false.
            do idom_loop = 1,dom_recv%size()
                if ( idomain_g == self%dom(idom_loop)%vecs(1)%dparent_g() ) then
                    dom_store = idom_loop
                    call self%dom_store%push_back(idom_loop)
                    domain_found = .true.
                    exit
                end if
            end do
            if ( .not. domain_found ) call chidg_signal(FATAL,"chidgVector_recv_comm%init: domain ordering, domain not found")





            do ielem_recv = 1,nelem_recv

                call MPI_Recv(ielement_g, 1, MPI_INTEGER4, self%proc, 0, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)

                ! Find element store index
                element_found = .false.
                do ielem_loop = 1,nelem_recv
                    if ( ielement_g == self%dom(dom_store)%vecs(ielem_loop)%eparent_g() ) then
                        call self%elem_store(idom_recv)%push_back(ielem_loop)
                        element_found = .true.
                        exit
                    end if
                end do
                if (.not. element_found) call chidg_signal(FATAL,"chidgVector_recv_comm%init: element ordering, element not found")

            end do ! ielem_recv


        end do ! idom_recv





    end subroutine init
    !*********************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/7/2016
    !!
    !!
    !----------------------------------------------------------------------------------
    subroutine clear(self)
        class(chidgVector_recv_comm_t), intent(inout)   :: self

        integer(ik) :: idom


        do idom = 1,size(self%dom)
            call self%dom(idom)%clear()
        end do


    end subroutine clear
    !**********************************************************************************







end module type_chidgVector_recv_comm
