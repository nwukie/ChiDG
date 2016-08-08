module type_RASILU0_send_comm
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use mod_constants,              only: NFACES, DIAG, INTERIOR
    use mod_chidg_mpi,              only: ChiDG_COMM, IRANK
    use type_ivector,               only: ivector_t
    use type_mesh,                  only: mesh_t
    use type_chidgMatrix,           only: chidgMatrix_t
    use type_RASILU0_send_comm_dom, only: RASILU0_send_comm_dom_t

    use mpi_f08,                    only: MPI_ISend, MPI_INTEGER4, MPI_REQUEST
    implicit none



    !>
    !!
    !!
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------------
    type, public :: RASILU0_send_comm_t

        integer(ik)                                 :: proc
        type(RASILU0_send_comm_dom_t),  allocatable :: dom(:)

    contains

        procedure   :: init

    end type RASILU0_send_comm_t
    !************************************************************************************************





contains



    !>  
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/22/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine init(self,mesh,A,proc)
        class(RASILU0_send_comm_t), intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh(:)
        type(chidgMatrix_t),        intent(in)      :: A
        integer(ik),                intent(in)      :: proc

        type(ivector_t)             :: dom_send
        integer(ik)                 :: idom, idom_send, ierr, iblk, ielem, iface, nelem_send, ielem_send, ielem_n, iface_n, iblk_send, nblks
        integer(ik),    allocatable :: send_procs_dom(:)
        logical                     :: comm_domain, overlap_elem
        type(MPI_REQUEST)           :: null_request

        !
        ! Set send processor
        !
        self%proc = proc



        !
        ! Accumulate the domains that send to processor 'proc'
        !
        do idom = 1,size(mesh)

            send_procs_dom = mesh(idom)%get_send_procs_local()

            ! Check if this domain is communicating with 'proc'
            comm_domain = any(send_procs_dom == proc)

            if (comm_domain) call dom_send%push_back(idom)

        end do !idom





        !
        ! Allocate storage for element indices being sent from each dom in dom_send
        !
        allocate(self%dom(dom_send%size()), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! For each domain sending info to 'proc', send data from chidgMatrix.
        !
        do idom_send = 1,dom_send%size()


            idom = dom_send%at(idom_send)
            self%dom(idom_send)%idomain_g = mesh(idom)%idomain_g
            self%dom(idom_send)%idomain_l = mesh(idom)%idomain_l

            !
            ! Loop through element faces and find neighbors that are off-processor on 'proc' to determine which elements to send as overlap data
            !
            do ielem = 1,mesh(idom)%nelem
                do iface = 1,NFACES

                    overlap_elem = (mesh(idom)%faces(ielem,iface)%ineighbor_proc == proc)

                    if (overlap_elem) then
                        call self%dom(idom_send)%elem_send%push_back_unique(ielem)
                        exit
                    end if

                end do !iface
            end do !ielem


            nelem_send = self%dom(idom_send)%elem_send%size()
            !allocate(self%dom(idom)%blk_send(nelem_send), stat=ierr)
            allocate(self%dom(idom_send)%blk_send(nelem_send), stat=ierr)
            if (ierr /= 0) call AllocationError


            !
            ! Communicate the number of elements being sent to 'proc' for idomain_g
            !
            call MPI_ISend(self%dom(idom_send)%elem_send%size_, 1, MPI_INTEGER4, proc, self%dom(idom_send)%idomain_g, ChiDG_COMM, null_request, ierr)





            !
            ! For each element in the overlap, determine which linearization blocks to send
            !
            do ielem_send = 1,self%dom(idom_send)%elem_send%size()
                ielem = self%dom(idom_send)%elem_send%at(ielem_send)

                !
                ! Search for blocks to send that couple with the off-processor domain
                !
                do iface = 1,NFACES

                    if ( (mesh(idom)%faces(ielem,iface)%ftype == INTERIOR)  .and. (mesh(idom)%faces(ielem,iface)%ineighbor_proc == proc) ) then
                        call self%dom(idom_send)%blk_send(ielem_send)%push_back(iface)
                    end if

                end do






                !
                ! Search neighbors to see if any of them are also overlapping blocks, because we would need to send their linearization as well.
                !
                do iface = 1,NFACES

                    ! Get neighbor for iface
                    if ( (mesh(idom)%faces(ielem,iface)%ftype == INTERIOR) .and. (mesh(idom)%faces(ielem,iface)%ineighbor_proc == IRANK) ) then
                        ielem_n = mesh(idom)%faces(ielem,iface)%ineighbor_element_l


                        ! Loop through the faces of the neighbor element to see if it is also in the RAS overlap. If so, add its index to the list of linearization blocks to send for ielem.
                        do iface_n = 1,NFACES
                            if (  (mesh(idom)%faces(ielem_n,iface_n)%ineighbor_proc == proc) .and. (mesh(idom)%faces(ielem_n,iface_n)%ineighbor_domain_g == self%dom(idom_send)%idomain_g)  ) then
                                call self%dom(idom_send)%blk_send(ielem_send)%push_back(iface)
                                exit
                            end if
                        end do

                    end if


                end do



                !
                ! Add the block diagonal to the list to send
                !
                call self%dom(idom_send)%blk_send(ielem_send)%push_back(DIAG)




                ! Communicate the number of blocks being sent to 'proc' for element
                call MPI_ISend(self%dom(idom_send)%blk_send(ielem_send)%size_, 1, MPI_INTEGER4, proc, self%dom(idom_send)%idomain_g, ChiDG_COMM, null_request, ierr)

                ! Communicate which blocks are being send to 'proc' for element
                nblks = self%dom(idom_send)%blk_send(ielem_send)%size()
                call MPI_ISend(self%dom(idom_send)%blk_send(ielem_send)%data_(1:nblks), nblks, MPI_INTEGER4, proc, self%dom(idom_send)%idomain_g, ChiDG_COMM, null_request, ierr)

                ! For each block send block initialization data
                do iblk_send = 1,self%dom(idom_send)%blk_send(ielem_send)%size()
                    iblk = self%dom(idom_send)%blk_send(ielem_send)%at(iblk_send)

                    call MPI_ISend(A%dom(idom)%lblks(ielem,iblk)%nrows_,       1, MPI_INTEGER4, proc, self%dom(idom_send)%idomain_g, ChiDG_COMM, null_request, ierr)
                    call MPI_ISend(A%dom(idom)%lblks(ielem,iblk)%ncols_,       1, MPI_INTEGER4, proc, self%dom(idom_send)%idomain_g, ChiDG_COMM, null_request, ierr)
                    call MPI_ISend(A%dom(idom)%lblks(ielem,iblk)%dparent_g_,   1, MPI_INTEGER4, proc, self%dom(idom_send)%idomain_g, ChiDG_COMM, null_request, ierr)
                    call MPI_ISend(A%dom(idom)%lblks(ielem,iblk)%dparent_l_,   1, MPI_INTEGER4, proc, self%dom(idom_send)%idomain_g, ChiDG_COMM, null_request, ierr)
                    call MPI_ISend(A%dom(idom)%lblks(ielem,iblk)%eparent_g_,   1, MPI_INTEGER4, proc, self%dom(idom_send)%idomain_g, ChiDG_COMM, null_request, ierr)
                    call MPI_ISend(A%dom(idom)%lblks(ielem,iblk)%eparent_l_,   1, MPI_INTEGER4, proc, self%dom(idom_send)%idomain_g, ChiDG_COMM, null_request, ierr)
                    call MPI_ISend(A%dom(idom)%lblks(ielem,iblk)%parent_proc_, 1, MPI_INTEGER4, proc, self%dom(idom_send)%idomain_g, ChiDG_COMM, null_request, ierr)
                end do !iblk_send


            end do !ielem send





        end do !idom send





    end subroutine init
    !*************************************************************************************************






end module type_RASILU0_send_comm
