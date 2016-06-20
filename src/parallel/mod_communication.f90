module mod_communication
    use mod_kinds,      only: rk, ik
    use type_mesh,      only: mesh_t
    use mod_chidg_mpi,  only: IRANK, NRANK
    use mpi_f08
    implicit none

!#include "mpif.h"







contains




    !>  Establish element/face neighbor communication.
    !!
    !!  For a given mesh instance for a domain, attempt to find neighboring elements with which
    !!  to communicate. This searches both in the same domain on the current processor(local), and also
    !!  across processors, where part of a domain may have been split onto another processor(global)
    !!
    !!  This does NOT establish Chimera communication.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/16/2016
    !!
    !!  @param[inout]   mesh    An array of mesh types on the local processor.
    !!
    !-------------------------------------------------------------------------------
    subroutine establish_communication(mesh)
        type(mesh_t),   intent(inout)   :: mesh(:)

        integer(ik) :: imesh, iproc, idomain_g, idomain_l, ielem_l, ierr
        integer(ik) :: ineighbor_domain_g, ineighbor_domain_l, ineighbor_element_g, ineighbor_element_l
        integer(ik) :: data(4), corner_indices(4)
        logical     :: has_domain, searching, neighbor_element
        logical     :: includes_corner_one, includes_corner_two, includes_corner_three, includes_corner_four


        !
        ! All ranks initialize local communication
        !
        do imesh = 1,size(mesh)

            call mesh(imesh)%init_comm_local()
            print*, "Rank ", IRANK, ": Established local communication"

        end do











        !
        ! Loop through mesh types
        !
        do iproc = 0,NRANK-1



            ! For current rank, send out requests for neighbors
            if ( iproc == IRANK ) then
                do imesh = 1,size(mesh)

                    print*, "Rank ", IRANK, ": Establishing global communication"
                    call mesh(imesh)%init_comm_global()
                    print*, "Rank ", IRANK, ": Done Establishing global communication"

                end do !imesh
            end if ! iproc




            ! 
            ! All other procs, listen for requests and send back matches
            !
            if ( iproc /= IRANK ) then

                ! Process face request and send status
                print*, "Rank ", IRANK, ": Waiting for request from root process"
                call MPI_BCast(searching,1,MPI_LOGICAL,iproc,MPI_COMM_WORLD,ierr)
                print*, "Rank ", IRANK, ": Done waiting for request from root process"
                do while ( searching )
                    
                    ! Receive global domain index of mesh being searched
                    print*, "Rank ", IRANK, ": Waiting for domain index"
                    call MPI_BCast(idomain_g,1,MPI_INTEGER4,iproc,MPI_COMM_WORLD,ierr)
                    print*, "Rank ", IRANK, ": Done waiting for domain index"

                    ! Check if we have domain index
                    has_domain = .false.
                    do imesh = 1,size(mesh)
                        if ( idomain_g == mesh(imesh)%idomain_g ) then
                            idomain_l  = imesh
                            has_domain = .true.
                        end if
                    end do

                    ! Send domain status
                    print*, "Rank ", IRANK, ": Sending domain status"
                    call MPI_Send(has_domain,1,MPI_LOGICAL,iproc,1,MPI_COMM_WORLD,ierr)
                    print*, "Rank ", IRANK, ": Done sending domain status"


                    ! If we have partition of the global domain being searched, receive face information and try to find match
                    if ( has_domain ) then 
                        ! Receive corner indices of face to be matched
                        print*, "Rank ", IRANK, ": Waiting for face corner indices"
                        call MPI_Recv(corner_indices,4,MPI_INTEGER4,iproc,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
                        print*, "Rank ", IRANK, ": Done waiting for face corner indices"

                        ! Loop through local domain and try to find a match
                        ! Test the incoming face nodes against local elements, if all face nodes are also contained in an element, then they are neighbors.
                        neighbor_element = .false.
                        do ielem_l = 1,mesh(idomain_l)%nelem
                            includes_corner_one   = any( mesh(idomain_l)%elems(ielem_l)%connectivity%get_element_nodes() == corner_indices(1) )
                            includes_corner_two   = any( mesh(idomain_l)%elems(ielem_l)%connectivity%get_element_nodes() == corner_indices(2) )
                            includes_corner_three = any( mesh(idomain_l)%elems(ielem_l)%connectivity%get_element_nodes() == corner_indices(3) )
                            includes_corner_four  = any( mesh(idomain_l)%elems(ielem_l)%connectivity%get_element_nodes() == corner_indices(4) )
                            neighbor_element = ( includes_corner_one .and. includes_corner_two .and. includes_corner_three .and. includes_corner_four )

                            ! Send element-found status. If found, send element index information.
                            print*, "Rank ", IRANK, ": Sending face status"
                            call MPI_Send(neighbor_element,1,MPI_LOGICAL,iproc,3,MPI_COMM_WORLD,ierr)
                            print*, "Rank ", IRANK, ": Done sending face status"
                            if ( neighbor_element ) then
                                ineighbor_domain_g  = mesh(idomain_l)%elems(ielem_l)%connectivity%get_domain_index()
                                ineighbor_domain_l  = idomain_l
                                ineighbor_element_g = mesh(idomain_l)%elems(ielem_l)%connectivity%get_element_index()
                                ineighbor_element_l = ielem_l

                                data = [ineighbor_domain_g, ineighbor_domain_l, ineighbor_element_g, ineighbor_element_l]
                                print*, "Rank ", IRANK, ": Sending face information"
                                call MPI_Send(data,4,MPI_INTEGER4,iproc,4,MPI_COMM_WORLD,ierr)
                                print*, "Rank ", IRANK, ": Done sending face information"
                                exit
                            end if

                        end do

                    end if

                    ! Detect if iproc is still sending out face requests. 
                    print*, "Rank ", IRANK, ": Waiting for request status from root process"
                    call MPI_BCast(searching,1,MPI_LOGICAL,iproc,MPI_COMM_WORLD,ierr)
                    print*, "Rank ", IRANK, ": Done waiting for request status from root process"
                end do 
                

            end if !iproc



            !
            ! Barrier
            !
            print*, "Rank ", IRANK, ": Waiting at Barrier"
            call MPI_Barrier(MPI_COMM_WORLD,ierr)
            print*, "Rank ", IRANK, ": Done waiting at Barrier"


        end do






    end subroutine establish_communication
    !********************************************************************************






























end module mod_communication
