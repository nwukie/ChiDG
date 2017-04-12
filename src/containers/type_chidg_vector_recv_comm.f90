module type_chidg_vector_recv_comm
#include <messenger.h>
    use mod_kinds,          only: ik
    use mod_constants,      only: INTERIOR, CHIMERA
    use type_mesh,          only: mesh_t
    use type_ivector,       only: ivector_t
    use type_domain_vector, only: domain_vector_t
    use mod_chidg_mpi,      only: ChiDG_COMM
    use mpi_f08,            only: MPI_Recv, MPI_INTEGER4, MPI_STATUS_IGNORE, MPI_ANY_TAG
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !!
    !----------------------------------------------------------------------------
    type, public :: chidg_vector_recv_comm_t

        integer(ik)                         :: proc
        type(domain_vector_t),  allocatable :: dom(:)

    contains

        procedure,  public  :: init
        procedure,  public  :: clear

    end type chidg_vector_recv_comm_t
    !****************************************************************************








contains




    !>  Receive information from 'proc' about what it is sending, so we can prepare
    !!  here about what to receive.
    !!
    !!  Communication here is initiated from 'type_chidg_vector_send_comm'.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/30/2016
    !!
    !!
    !---------------------------------------------------------------------------------
    subroutine init(self,mesh,proc,comm)
        class(chidg_vector_recv_comm_t), intent(inout)   :: self
        type(mesh_t),                    intent(inout)   :: mesh
        integer(ik),                     intent(in)      :: proc
        integer(ik),                     intent(in)      :: comm

        character(:),   allocatable :: user_msg
        integer(ik)                 :: idom, idom_recv, ndom_recv, ierr, ielem, iface,  &
                                       ChiID, donor_domain_g, donor_element_g,          &
                                       idomain_g, ielement_g, dom_store, idom_loop,     &
                                       ielem_loop, ielem_recv, nelem_recv,              &
                                       neighbor_domain_g, neighbor_element_g,           &
                                       recv_element, recv_domain, idonor
        integer(ik),    allocatable :: comm_procs_dom(:)
        logical                     :: proc_has_domain, domain_found, element_found,    &
                                       comm_neighbor, is_interior, is_chimera,         &
                                       comm_donor, donor_recv_found

        !
        ! Set processor being received from
        !
        self%proc = proc


        
        !
        ! Get number of domains being received from proc. This is sent from chidg_vector_send_comm
        !
        call MPI_Recv(ndom_recv,1,MPI_INTEGER4,self%proc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)



        !
        ! Allocate number of recv domains to recv information from proc
        !
        if (allocated(self%dom)) deallocate(self%dom)
        allocate(self%dom(ndom_recv), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Call initialization for each domain to be received
        !
        do idom_recv = 1,ndom_recv
            call self%dom(idom_recv)%init(proc)
        end do ! idom





        !
        ! Loop through mesh and initialize recv indices
        !
        do idom = 1,mesh%ndomains()
            do ielem = 1,mesh%domain(idom)%nelem
                do iface = 1,size(mesh%domain(idom)%faces,2)

                    is_interior = ( mesh%domain(idom)%faces(ielem,iface)%ftype == INTERIOR )
                    is_chimera  = ( mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA  )


                    !
                    ! Initialize recv for comm interior neighbors
                    !
                    if (is_interior) then

                        comm_neighbor = (proc == mesh%domain(idom)%faces(ielem,iface)%ineighbor_proc)

                        ! If neighbor is being communicated, find it in the recv domains
                        if (comm_neighbor) then
                            neighbor_domain_g  = mesh%domain(idom)%faces(ielem,iface)%ineighbor_domain_g
                            neighbor_element_g = mesh%domain(idom)%faces(ielem,iface)%ineighbor_element_g

                            ! Loop through domains being received to find the right domain
                            do idom_recv = 1,ndom_recv
                                recv_domain = self%dom(idom_recv)%vecs(1)%dparent_g()

                                if (recv_domain == neighbor_domain_g) then

                                    ! Loop through the elements in the recv domain to find the right neighbor element
                                    do ielem_recv = 1,size(self%dom(idom_recv)%vecs)
                                        recv_element = self%dom(idom_recv)%vecs(ielem_recv)%eparent_g()

                                        

                                        ! Set the location where a face can find its off-processor neighbor 
                                        if (recv_element == neighbor_element_g) then
                                            mesh%domain(idom)%faces(ielem,iface)%recv_comm    = comm
                                            mesh%domain(idom)%faces(ielem,iface)%recv_domain  = idom_recv
                                            mesh%domain(idom)%faces(ielem,iface)%recv_element = ielem_recv
                                        end if

                                    end do !ielem_recv

                                end if

                            end do !idom_recv

                        end if





                    !
                    ! Initialize recv for comm chimera receivers
                    !
                    else if (is_chimera) then
                        
                        ChiID = mesh%domain(idom)%faces(ielem,iface)%ChiID
                        do idonor = 1,mesh%domain(idom)%chimera%recv%data(ChiID)%ndonors()

                            comm_donor = (proc == mesh%domain(idom)%chimera%recv%data(ChiID)%donor_proc%at(idonor) )

                            donor_recv_found = .false.
                            if (comm_donor) then
                                donor_domain_g  = mesh%domain(idom)%chimera%recv%data(ChiID)%donor_domain_g%at(idonor)
                                donor_element_g = mesh%domain(idom)%chimera%recv%data(ChiID)%donor_element_g%at(idonor)


                                ! Loop through domains being received to find the right domain
                                do idom_recv = 1,ndom_recv
                                    recv_domain = self%dom(idom_recv)%vecs(1)%dparent_g()

                                    if (recv_domain == donor_domain_g) then

                                        ! Loop through the elements in the recv domain to find the right neighbor element
                                        do ielem_recv = 1,size(self%dom(idom_recv)%vecs)
                                            recv_element = self%dom(idom_recv)%vecs(ielem_recv)%eparent_g()


                                            ! Set the location where a face can find its off-processor neighbor 
                                            if (recv_element == donor_element_g) then
                                                mesh%domain(idom)%chimera%recv%data(ChiID)%donor_recv_comm%data_(idonor)    = comm
                                                mesh%domain(idom)%chimera%recv%data(ChiID)%donor_recv_domain%data_(idonor)  = idom_recv
                                                mesh%domain(idom)%chimera%recv%data(ChiID)%donor_recv_element%data_(idonor) = ielem_recv
                                                donor_recv_found = .true.
                                                exit
                                            end if

                                        end do !ielem_recv

                                    end if

                                    if (donor_recv_found) exit
                                end do !idom_recv

                                user_msg = "chidg_vector_recv_comm%init: chimera receiver did not &
                                            find parallel donor."
                                if (.not. donor_recv_found) call chidg_signal_three(FATAL,user_msg,IRANK,donor_domain_g,donor_element_g)
                            end if !comm_donor


                        end do ! idonor

                    

                    end if



                end do ! iface
            end do ! ielem
        end do ! idom





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
        class(chidg_vector_recv_comm_t), intent(inout)   :: self

        integer(ik) :: idom


        do idom = 1,size(self%dom)
            call self%dom(idom)%clear()
        end do


    end subroutine clear
    !**********************************************************************************







end module type_chidg_vector_recv_comm
