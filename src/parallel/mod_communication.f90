module mod_communication
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use type_mesh,  only: mesh_t
    use mod_chidg_mpi,  only: IRANK, NRANK, GLOBAL_MASTER
    use mod_chimera,    only: detect_chimera_faces,  &
                              detect_chimera_donors, &
                              compute_chimera_interpolators
    use mpi_f08
    implicit none

contains




    !>  Establish element/face neighbor communication.
    !!
    !!  For a given mesh instance for a domain, attempt to find neighboring elements with which
    !!  to communicate. This searches both in the same domain on the current processor(local), 
    !!  and also across processors, where part of a domain may have been split onto another 
    !!  processor(global)
    !!
    !!  This does NOT establish Chimera communication.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/16/2016
    !!
    !!  @param[inout]   mesh        An array of mesh types on the local processor.
    !!  @param[in]      ChiDG_COMM  An mpi communicator of the relevant processors. 
    !!                              Particularly usefull for testing.
    !!
    !------------------------------------------------------------------------------------------
    subroutine establish_neighbor_communication(mesh,ChiDG_COMM)
        type(mesh_t),   intent(inout)   :: mesh
        type(mpi_comm),     intent(in)      :: ChiDG_COMM

        integer(ik) :: idom, iproc, idomain_g, idomain_l, ielem_l, ierr
        integer(ik) :: ineighbor_domain_g, ineighbor_domain_l, &
                       ineighbor_element_g, ineighbor_element_l
        integer(ik) :: data(4), corner_indices(4)
        logical     :: has_domain, searching, neighbor_element, searching_mesh
        logical     :: includes_corner_one, includes_corner_two, &
                       includes_corner_three, includes_corner_four


        !
        ! All ranks initialize local communication
        !
        call write_line("Initialize: proc-local neighbor communication...", io_proc=GLOBAL_MASTER)
        do idom = 1,mesh%ndomains()
            call mesh%domain(idom)%init_comm_local()
        end do






        !
        ! Loop through mesh types
        !
        call write_line("Initialize: proc-global neighbor communication...", io_proc=GLOBAL_MASTER)
        do iproc = 0,NRANK-1



            !
            ! For current rank, send out requests for neighbors
            !
            if ( iproc == IRANK ) then

                do idom = 1,mesh%ndomains()
                    ! Broadcast that a mesh from iproc is searching for face neighbors
                    searching_mesh=.true.
                    call MPI_BCast(searching_mesh,1,MPI_LOGICAL,iproc,ChiDG_COMM,ierr)

                    ! Initiate search for communication across processors
                    call mesh%domain(idom)%init_comm_global(ChiDG_COMM)

                end do !idom


                ! Broadcast that no more mesh instances are searching for face neighbors
                searching_mesh=.false.
                call MPI_BCast(searching_mesh,1,MPI_LOGICAL,iproc,ChiDG_COMM,ierr)
            end if ! iproc





            ! 
            ! All other procs, listen for requests and send back matches
            !
            if ( iproc /= IRANK ) then


                ! Check if iproc is searching mesh instances
                call MPI_BCast(searching_mesh,1,MPI_LOGICAL,iproc,ChiDG_COMM,ierr)
                do while ( searching_mesh )


                    call MPI_BCast(searching,1,MPI_LOGICAL,iproc,ChiDG_COMM,ierr)
                    do while ( searching )
                        
                        ! Receive global domain index of mesh being searched
                        call MPI_Recv(idomain_g,1,MPI_INTEGER4,iproc,0,ChiDG_COMM,MPI_STATUS_IGNORE,ierr)


                        ! Search through local domains and check indices
                        do idom = 1,mesh%ndomains()
                            has_domain = ( idomain_g == mesh%domain(idom)%idomain_g )
                            if ( has_domain ) then
                                exit
                            end if
                        end do


                        
                        ! Send domain status. If we have a match, initiate handle_neighbor_request
                        call MPI_Send(has_domain,1,MPI_LOGICAL,iproc,1,ChiDG_COMM,ierr)

                        if ( has_domain ) then
                            call mesh%domain(idom)%handle_neighbor_request(iproc,ChiDG_COMM)
                        end if



                        ! Detect if iproc is still sending out face requests. 
                        call MPI_BCast(searching,1,MPI_LOGICAL,iproc,ChiDG_COMM,ierr)

                    end do ! searching faces


                    ! Check if iproc is searching another mesh or finished.
                    call MPI_BCast(searching_mesh,1,MPI_LOGICAL,iproc,ChiDG_COMM,ierr)
                end do ! searching_mesh

            end if !iproc





            !
            ! Barrier
            !
            call MPI_Barrier(ChiDG_COMM,ierr)


        end do ! iproc






    end subroutine establish_neighbor_communication
    !*****************************************************************************************



















    !>  Establish chimera face communication.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/22/2016
    !!
    !!
    !!  @param[inout]   mesh        An array of mesh types on the local processor.
    !!  @param[in]      ChiDG_COMM  An mpi communicator of the relevant processors. 
    !!                              Particularly useful for testing.
    !!
    !------------------------------------------------------------------------------------------
    subroutine establish_chimera_communication(mesh, ChiDG_COMM)
        type(mesh_t),   intent(inout)   :: mesh
        type(mpi_comm)                      :: ChiDG_COMM

        integer :: idom, ierr


        !
        ! In case establish_chimera_communication is being called as a reinitialization 
        ! procedure. Calling chimera%clear to wipe out previous data 
        ! and redetect all chimera faces, donors and reinitialize donor data.
        !
        call write_line("Initialize: chimera communication...", io_proc=GLOBAL_MASTER)
        do idom = 1,mesh%ndomains()
            !call mesh%domain(idom)%chimera%send%clear()
            !call mesh%domain(idom)%chimera%recv%clear()
            call mesh%domain(idom)%chimera%clear()
        end do !idom
        call MPI_Barrier(ChiDG_COMM,ierr)   ! not sure if this is needed.



        ! Process-local routine, just to flag faces
        call detect_chimera_faces(mesh)


        ! Detect local or global chimera donors
        call detect_chimera_donors(mesh)


        ! Compute chimera interpolation matrices to evaluate donor 
        ! solution at receiver quadrature nodes.
        call compute_chimera_interpolators(mesh)


        ! Barrier
        call MPI_Barrier(ChiDG_COMM,ierr)

    end subroutine establish_chimera_communication
    !*****************************************************************************************










end module mod_communication
