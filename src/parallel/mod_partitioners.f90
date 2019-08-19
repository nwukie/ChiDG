module mod_partitioners
#include <messenger.h>
    use type_meshdata,              only: meshdata_t
    use mod_chidg_mpi,              only: NRANK, IRANK, GLOBAL_MASTER
    use mod_metis,                  only: METIS_PartMeshNodal, METIS_PartMeshDual
    use type_mpi_request_vector,    only: mpi_request_vector_t
    use mpi_f08

    use type_domain_connectivity,   only: domain_connectivity_t
    use type_partition,             only: partition_t
    implicit none



    type(mpi_request_vector_t)  :: send_partition_requests



contains




    !>  Take connectivity information from multiple domains, combine them into one
    !!  global connectivity, and partition.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/8/2016
    !! 
    !!
    !!  @param[in]      connectivities  Array of domain connectivity types, one for each domain
    !!  @param[inout]   partitions      Mesh partition data returned
    !!
    !-----------------------------------------------------------------------------------------
    subroutine partition_connectivity(connectivities, weights, partitions)
        use iso_c_binding,  only: c_ptr, c_int, c_null_ptr
        type(domain_connectivity_t),    intent(inout)   :: connectivities(:)
        real(rk),                       intent(in)      :: weights(:)
        type(partition_t), allocatable, intent(inout)   :: partitions(:)

        character(:),   allocatable :: domain_name
        logical                     :: serial, parallel
        logical                     :: element_in_partition
        integer(ik),    allocatable :: npartition_elements_in_domain(:), element_nodes(:)
        integer(ik)                 :: ielem_part, ipart_conn, ndomains_in_partition, ipart_elem

        integer                     :: idomain, ndomains, ielem_conn, nelem_conn, nelem_total, &
                                       ierr, idomain_g, ielement_g, ielem, inode, ielem_node,  &
                                       node_offset, iconn, nconn, mapping, nnodes_element,     &
                                       nnodes_conn, ipartition, nindices_total,                &
                                       nindices_element, element_partition

        integer(c_int)                      :: nelem, nnodes
        integer(c_int), allocatable         :: eptr(:), eind(:), epart(:), npart(:)
        integer(c_int), allocatable, target :: vwgt(:), options(:)
        type(c_ptr)                         :: vwgt_p, options_p, vsize, tpwgts
        integer(c_int)                      :: n, npartitions, ncommon


        serial   = (NRANK == 1)
        parallel = (NRANK > 1)

        !
        ! Allocate partitions
        !
        allocate(partitions(NRANK), stat=ierr)
        if (ierr /= 0) call AllocationError



        if ( serial ) then

            ! All domains on one partition
            nconn = size(connectivities)
            call partitions(1)%init(nconn)
            partitions(1)%connectivities = connectivities        

            ! Set default partition
            do iconn = 1,size(partitions(1)%connectivities)
                nelem = partitions(1)%connectivities(iconn)%get_nelements()
                do ielem = 1,nelem
                    call partitions(1)%connectivities(iconn)%data(ielem)%set_element_partition(IRANK)
                end do
            end do

        elseif ( parallel ) then


            ! Unused METIS parameters
            vwgt_p  = c_null_ptr
            vsize   = c_null_ptr
            tpwgts  = c_null_ptr

            ! METIS Options
            allocate(options(0:40),stat=ierr)
            if (ierr /= 0) call AllocationError

            call METIS_SetDefaultOptions(options)
            options_p = c_loc(options(0))

            options(0)  = 0      ! PTYPE  (Recursive Bisetion (a), k-way (b)
            options(3)  = 0      ! IPTYPE (

            !options(10) = 2
            !options(16) = 1

            !
            ! Get number of domains
            !
            nconn = size(connectivities)


            !
            ! Accumulate total number of elements and nodes across all domains
            !
            nelem_total    = 0
            nindices_total = 0
            nnodes         = 0
            do iconn = 1,nconn
                ! Accumulate total number of elements across domains
                nelem_conn = connectivities(iconn)%get_nelements()
                nelem_total  = nelem_total  + nelem_conn


                ! Accumulate total number of element indices across domains
                do ielem_conn = 1,nelem_conn
                    mapping          = connectivities(iconn)%get_element_mapping(ielem_conn)
                    nindices_element = (mapping+1)*(mapping+1)*(mapping+1)
                    nindices_total   = nindices_total + nindices_element
                end do

                nnodes = nnodes + connectivities(iconn)%get_nnodes()
            end do


            !
            ! Allocate storage for partitioning
            !
            allocate(eptr(nelem_total+1), eind(nindices_total), vwgt(nelem_total), epart(nelem_total), npart(nnodes), stat=ierr)
            if (ierr /= 0) call AllocationError





            !
            ! Assemble connectivities into METIS-readable format
            !
            inode = 1
            ielem = 1
            node_offset = 0
            do iconn = 1,nconn

                nelem_conn = connectivities(iconn)%get_nelements()

                do ielem_conn = 1,nelem_conn

                    mapping = connectivities(iconn)%get_element_mapping(ielem_conn)
                    nnodes_element = (mapping+1)*(mapping+1)*(mapping+1)
                    ! Set pointer index to beginning of element connectivity
                    if ( iconn == 1 .and. ielem_conn == 1) then
                        eptr(ielem) = 0
                    else
                        eptr(ielem) = eptr(ielem-1) + nnodes_element
                    end if


                    ! Set element connectivity
                    do ielem_node = 1,nnodes_element
                        eind(inode) = connectivities(iconn)%get_element_node(ielem_conn,ielem_node) + node_offset - 1    ! 0-based
                        inode = inode + 1
                    end do

                    ! Set element weight
                    !vwgt(ielem) = domain_weights(iconn)
                    vwgt(ielem) = 1

                    ielem = ielem + 1
                end do ! ielem_domain
            

                ! Offset by number of nodes in all previous domains
                node_offset = node_offset + connectivities(iconn)%get_nnodes()

            end do !iconn





            !
            ! METIS documentation is not clear why/what this last value is for.
            !
            eptr(nelem_total+1) = size(eind)








            !
            ! Partition mesh using METIS
            !
            npartitions = NRANK
            ncommon     = int(4,c_int)

            !call METIS_PartMeshNodal(nelem_total,nnodes,eptr,eind,vwgt,vsize,npartitions,tpwgts,options,n,epart,npart)
            !ierr = METIS_PartMeshNodal(nelem_total,nnodes,eptr,eind,vwgt,vsize,npartitions,tpwgts,options,n,epart,npart)

            ! Don't seem to be able to get the node weights correct with this so higher-order elements are weighted appropriately.
            ! PartMeshDual seems to do the trick. Seems to partition by elements instead of nodes.
            !ierr = METIS_PartMeshNodal(nelem_total,nnodes,eptr,eind,vwgt,vsize,npartitions,tpwgts,options_p,n,epart,npart)


            ! Standard
            vwgt_p = c_loc(vwgt(1))
            ierr = METIS_PartMeshDual(nelem_total,nnodes,eptr,eind,vwgt_p,vsize,ncommon,npartitions,tpwgts,options_p,n,epart,npart)
            if (ierr /= 1) call chidg_signal(FATAL,"partition_connectivity: Error in METIS_PartMeshNodal")





            !
            ! Distribute global partition data to domain-level connectivities
            !
            ielem = 1
            do iconn = 1,nconn
                nelem_conn = connectivities(iconn)%get_nelements()


                do ielem_conn = 1,nelem_conn
                    ! Distribute partition indices from global array to domain-level connectivities
                    call connectivities(iconn)%data(ielem_conn)%set_element_partition(epart(ielem))

                    ielem = ielem + 1   ! global index
                end do ! ielem_domain


            end do ! iconn
            



            !
            ! For each partition, search for elements that belong in the connectivity arrays.
            !
            allocate(npartition_elements_in_domain(nconn), stat=ierr)
            if (ierr /= 0) call AllocationError



            do ipartition = 1,npartitions

                npartition_elements_in_domain = 0
                ndomains_in_partition         = 0
                




                !
                ! Loop through connectivities and detect number of elements(if any) belong to the current partitions
                !
                do iconn = 1,nconn
                    nelem_conn = connectivities(iconn)%get_nelements()
                        
                    do ielem_conn = 1,nelem_conn
                        element_partition = connectivities(iconn)%data(ielem_conn)%get_element_partition()
                        if ( element_partition  ==  ipartition-1 ) then
                            npartition_elements_in_domain(iconn) = npartition_elements_in_domain(iconn) + 1
                        end if
                    end do ! ielem_conn


                    if ( npartition_elements_in_domain(iconn) > 0 ) then
                        ndomains_in_partition = ndomains_in_partition + 1
                    end if

                end do !iconn




                !
                ! Initialize partition with ndomains
                !
                !if (ndomains_in_partition == 0) call chidg_signal(FATAL,"partition_connectivity: No partition generated by metis for one of the processes.")
                call partitions(ipartition)%init(ndomains_in_partition)

                


                !
                ! Loop through and collect connectivities from included domains
                !
                ipart_conn = 0
                do iconn = 1,nconn

                    ! Get domain index of the current connectivity
                    idomain = connectivities(iconn)%get_domain_index()


                    if ( npartition_elements_in_domain(iconn) > 0 ) then
                        ipart_conn = ipart_conn + 1


                        ! Initialize partition connectivity
                        domain_name = connectivities(iconn)%get_domain_name()
                        nelem = npartition_elements_in_domain(iconn)
                        nnodes = connectivities(iconn)%get_nnodes()
                        call partitions(ipartition)%connectivities(ipart_conn)%init(domain_name,nelem,nnodes)
                         
                        
                        ! Collect elements
                        nelem_conn = connectivities(iconn)%get_nelements()
                        ielem_part   = 1
                        do ielem_conn = 1,nelem_conn
                            element_in_partition = connectivities(iconn)%get_element_partition(ielem_conn) == (ipartition - 1)
                            if ( element_in_partition ) then
                                ! Get information from main connectivity
                                idomain_g       = connectivities(iconn)%data(ielem_conn)%get_domain_index()
                                ielement_g      = connectivities(iconn)%data(ielem_conn)%get_element_index()
                                mapping         = connectivities(iconn)%data(ielem_conn)%get_element_mapping()
                                ipart_elem      = connectivities(iconn)%data(ielem_conn)%get_element_partition()
                                element_nodes   = connectivities(iconn)%data(ielem_conn)%get_element_nodes()


                                ! Set information on element in partition
                                call partitions(ipartition)%connectivities(ipart_conn)%data(ielem_part)%init(mapping)
                                call partitions(ipartition)%connectivities(ipart_conn)%data(ielem_part)%set_domain_index(idomain_g)
                                call partitions(ipartition)%connectivities(ipart_conn)%data(ielem_part)%set_element_index(ielement_g)
                                call partitions(ipartition)%connectivities(ipart_conn)%data(ielem_part)%set_element_mapping(mapping)
                                call partitions(ipartition)%connectivities(ipart_conn)%data(ielem_part)%set_element_nodes(element_nodes)
                                call partitions(ipartition)%connectivities(ipart_conn)%data(ielem_part)%set_element_partition(ipart_elem)

                                ielem_part = ielem_part + 1
                            end if
                        end do

                    end if

                end do !iconn


            end do ! ipartition


        end if !serial/parallel





    end subroutine partition_connectivity
    !*************************************************************************************************














    !>  Master process sends partition information.
    !!
    !!  NOTE: take care with non-blocking sends, as the call assumes that the data in the 
    !!        location will remain valid at a future time. For example, can't send local 
    !!        variables in the subroutine since their memory is released once the subroutine 
    !!        ends. For this reason, all the sends are to variables within the paritions(:) 
    !!        containers. This persists.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine send_partitions(partitions,ChiDG_COMM)
        type(partition_t),  intent(inout), asynchronous  :: partitions(:)
        type(mpi_comm),     intent(in)  :: ChiDG_COMM

        integer                     :: ipartition, npartitions, ierr
        integer                     :: ielem, data_size
        integer(ik)                 :: iconn , nelements
        integer(ik)                 :: nconn
        type(MPI_Request)           :: handle, handle1, handle2, handle3, handle4


        npartitions = size(partitions)


        !
        ! Send each partition to its processor
        !
        do ipartition = 1,npartitions


            ! Send number of connectivities in current partition
            call MPI_ISend(partitions(ipartition)%nconn,1,MPI_INTEGER4, ipartition-1, 0, ChiDG_COMM, handle, ierr)
            call send_partition_requests%push_back(handle)
            


            ! Send connectivities
            nconn = partitions(ipartition)%nconn
            do iconn = 1,nconn

                !string_length = partitions(ipartition)%connectivities(iconn)%name_length
                !call MPI_ISend(partitions(ipartition)%connectivities(iconn)%name_length,  1,                    MPI_INTEGER4,  ipartition-1, 0, ChiDG_COMM, handle1, ierr)
                !call MPI_ISend(partitions(ipartition)%connectivities(iconn)%name,         int(string_length,4), MPI_CHARACTER, ipartition-1, 0, ChiDG_COMM, handle2, ierr)
                call MPI_ISend(partitions(ipartition)%connectivities(iconn)%name,         1000,                 MPI_CHARACTER, ipartition-1, 0, ChiDG_COMM, handle2, ierr)
                call MPI_ISend(partitions(ipartition)%connectivities(iconn)%nelements,    1,                    MPI_INTEGER4,  ipartition-1, 0, ChiDG_COMM, handle3, ierr)
                call MPI_ISend(partitions(ipartition)%connectivities(iconn)%nnodes,       1,                    MPI_INTEGER4,  ipartition-1, 0, ChiDG_COMM, handle4, ierr)


                !call send_partition_requests%push_back(handle1)
                call send_partition_requests%push_back(handle2)
                call send_partition_requests%push_back(handle3)
                call send_partition_requests%push_back(handle4)

                ! Send a connectivity for each element. Would be more efficient to assemble a consolodated array, but must
                ! be careful, because it can't be reused until it has been received since we are posting non-blocking sends.
                nelements = partitions(ipartition)%connectivities(iconn)%get_nelements()
                do ielem = 1,nelements
                    call MPI_ISend(partitions(ipartition)%connectivities(iconn)%data(ielem)%connectivity_size, 1, MPI_INTEGER4, ipartition-1, 0, ChiDG_COMM, handle1, ierr)
                    data_size = partitions(ipartition)%connectivities(iconn)%data(ielem)%connectivity_size
                    call MPI_ISend(partitions(ipartition)%connectivities(iconn)%data(ielem)%data, data_size, MPI_INTEGER4, ipartition-1, 0, ChiDG_COMM, handle2, ierr)

                    call send_partition_requests%push_back(handle1)
                    call send_partition_requests%push_back(handle2)
                end do

            end do ! iconn


        end do ! ipartition


    end subroutine send_partitions
    !*****************************************************************************************










    !>  Receives partition information from the master process.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/16/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine recv_partition(partition,ChiDG_COMM)
        type(partition_t),  intent(inout)   :: partition
        type(mpi_comm),     intent(in)      :: ChiDG_COMM

        integer                     :: ipartition, npartitions, ierr
        integer                     :: dims(3), iconn, idomain, ielem
        integer                     :: idomain_g, ielement_g, mapping, nnodes_element
        integer(ik), volatile                 :: conn_size, nnodes, nelements, nconn, nrequests
        integer,        allocatable :: nodes(:)
        integer(ik),    allocatable :: conn(:)
        character(1000)             :: domain_name
        character(:),   allocatable :: domain_name_trim



        ! Receive number of connectivities and inititalize partition
        call MPI_Recv(nconn,1,MPI_INTEGER4, GLOBAL_MASTER, 0, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
        call partition%init(nconn)


        ! Receive connectivities
        do iconn = 1,nconn

!            call MPI_Recv(string_length,1,MPI_INTEGER4, GLOBAL_MASTER, 0, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
!
!            if (allocated(domain_name)) deallocate(domain_name)
!            allocate(character(len=string_length) :: domain_name, stat=ierr)
!            if (ierr /= 0) call AllocationError

            !call MPI_Recv(domain_name,int(string_length,4),MPI_CHARACTER, GLOBAL_MASTER, 0, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
            call MPI_Recv(domain_name,1000,MPI_CHARACTER, GLOBAL_MASTER, 0, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
            domain_name_trim = trim(domain_name)



            ! Recv connectivity information
            call MPI_Recv(nelements,1,MPI_INTEGER4, GLOBAL_MASTER, 0, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
            call MPI_Recv(nnodes,   1,MPI_INTEGER4, GLOBAL_MASTER, 0, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)

            call partition%connectivities(iconn)%init(domain_name_trim,nelements,nnodes)



            ! Recv connectivity for each element
            do ielem = 1,nelements

                ! Get size to recv
                call MPI_Recv(conn_size, 1, MPI_INTEGER4, GLOBAL_MASTER, 0, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)

                if (allocated(conn)) deallocate(conn)
                allocate(conn(conn_size), stat=ierr)
                if (ierr /= 0) call AllocationError

                ! Recv connectivity
                call MPI_Recv(conn,conn_size, MPI_INTEGER4, GLOBAL_MASTER, 0, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)


                ! Set element connectivity
                idomain_g      = conn(1)
                ielement_g     = conn(2)
                mapping        = conn(3)
                nnodes_element = (mapping+1)*(mapping+1)*(mapping+1)
                nodes          = conn(4:4+nnodes_element-1)

                call partition%connectivities(iconn)%data(ielem)%init(mapping)
                call partition%connectivities(iconn)%data(ielem)%set_domain_index(idomain_g)
                call partition%connectivities(iconn)%data(ielem)%set_element_index(ielement_g)
                call partition%connectivities(iconn)%data(ielem)%set_element_mapping(mapping)
                call partition%connectivities(iconn)%data(ielem)%set_element_nodes(nodes)
                call partition%connectivities(iconn)%data(ielem)%set_element_partition(IRANK)

            end do !ielem


        end do !iconn




        !
        ! Wait on all requests initiated in send_partition
        !
        nrequests = send_partition_requests%size()
        if (nrequests > 0) then
            call MPI_Waitall(nrequests, send_partition_requests%data, MPI_STATUSES_IGNORE, ierr)
        end if



    end subroutine recv_partition
    !*****************************************************************************************




















end module mod_partitioners
