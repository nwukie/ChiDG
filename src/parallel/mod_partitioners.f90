module mod_partitioners
#include <messenger.h>
    use type_meshdata,  only: meshdata_t
    use mod_chidg_mpi,  only: NRANK, IRANK, GLOBAL_MASTER
    use mod_metis,      only: METIS_PartMeshNodal
    use mpi_f08

    use type_domain_connectivity,   only: domain_connectivity_t
    use type_partition,             only: partition_t
    implicit none







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
    !-----------------------------------------------------------------------------------------------
    subroutine partition_connectivity(connectivities, partitions)
        use iso_c_binding,  only: c_ptr, c_int, c_null_ptr
        type(domain_connectivity_t),    intent(inout) :: connectivities(:)
        type(partition_t), allocatable, intent(inout) :: partitions(:)

        logical                     :: serial, parallel

        logical                     :: element_in_partition
        integer(ik),    allocatable :: npartition_elements_in_domain(:), element_nodes(:)
        integer(ik)                 :: ielem_part, ipart_conn, ndomains_in_partition, ipart_elem

        integer                     :: idomain, ndomains, ielem_conn, nelem_conn, nelem_total, ierr, idomain_g, ielement_g
        integer                     :: ielem, inode, ielem_node, node_offset, iconn, nconn, mapping, nnodes_element
        integer                     :: ipartition, nindices_total, nindices_element, element_partition

        integer(c_int)              :: nelem, nnodes
        integer(c_int), allocatable :: eptr(:), eind(:), epart(:), npart(:)
        integer(c_int)              :: n, npartitions
        type(c_ptr)                 :: vwgt, vsize, options, tpwgts


        serial   = (NRANK == 1)
        parallel = (NRANK > 1)

        !
        ! Allocate partitions
        !
        allocate(partitions(NRANK), stat=ierr)
        if (ierr /= 0) call AllocationError



        if ( serial ) then
            ! All domains on one partition
            partitions(1)%connectivities = connectivities        


        elseif ( parallel ) then


            ! Unused METIS parameters
            vwgt    = c_null_ptr
            vsize   = c_null_ptr
            options = c_null_ptr
            tpwgts  = c_null_ptr


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
            allocate(eptr(nelem_total+1), eind(nindices_total), epart(nelem_total), npart(nnodes), stat=ierr)
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
            call METIS_PartMeshNodal(nelem_total,nnodes,eptr,eind,vwgt,vsize,npartitions,tpwgts,options,n,epart,npart)




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
                        nelem = npartition_elements_in_domain(iconn)
                        nnodes = connectivities(iconn)%get_nnodes()
                        call partitions(ipartition)%connectivities(ipart_conn)%init(nelem,nnodes)
                         
                        
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




        print*, 'number of partitions:', size(partitions)
        print*, '   number of connectivitites in partition 1:', size(partitions(1)%connectivities)
        !print*, '   number of connectivitites in partition 2:', size(partitions(2)%connectivities)



    end subroutine partition_connectivity
    !*************************************************************************************************














    !>  Master process sends partition information
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/8/2016
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine send_partitions(partitions)
        type(partition_t),  intent(in)    :: partitions(:)

        integer                     :: ipartition, npartitions, ierr
        integer                     :: dims(3), iconn, nconn, idomain, nelements, max_mapping, header_size, max_element_nodes, nnodes
        integer                     :: ielem
        integer(ik), allocatable    :: conn(:,:)


        npartitions = size(partitions)


        !
        ! Send each partition to its processor
        !
        do ipartition = 1,npartitions


            ! Send number of connectivities in current partition
            nconn = size(partitions(ipartition)%connectivities)
            print*, 'Number of connectivities in the partition being sent:', nconn
            call MPI_Send(nconn,1,MPI_INTEGER, ipartition-1, 0, MPI_COMM_WORLD, ierr)


            ! Send connectivities
            do iconn = 1,nconn

                ! Send connectivity information
                nelements         = partitions(ipartition)%connectivities(iconn)%get_nelements()
                max_mapping       = partitions(ipartition)%connectivities(iconn)%get_max_mapping()
                max_element_nodes = (max_mapping+1)*(max_mapping+1)*(max_mapping+1)
                header_size       = partitions(ipartition)%connectivities(iconn)%get_header_size()
                nnodes            = partitions(ipartition)%connectivities(iconn)%get_nnodes()
                dims(1) = nelements
                dims(2) = header_size + max_element_nodes
                dims(3) = nnodes
                call MPI_Send(dims,3,MPI_INTEGER, ipartition-1, 1+iconn, MPI_COMM_WORLD, ierr)

                ! Assemble connectivity data
                if (allocated(conn)) deallocate(conn)
                allocate(conn(dims(1),dims(2)), stat=ierr)
                if (ierr /= 0) call AllocationError
                do ielem = 1,nelements
                    conn(ielem,:) = partitions(ipartition)%connectivities(iconn)%data(ielem)%data
                end do

                ! Send connectivity data
                call MPI_Send(conn,dims(1)*dims(2), MPI_INTEGER4, ipartition-1, 1+nconn+1, MPI_COMM_WORLD, ierr)

            end do ! iconn



        end do ! ipartition


    end subroutine send_partitions
    !***************************************************************************************************










    !>  Receives partition information from the master process.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/16/2016
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine recv_partition(partition)
        type(partition_t),  intent(inout) :: partition

        integer                     :: ipartition, npartitions, ierr
        integer                     :: dims(3), iconn, nconn, idomain, nelements, ielem, iprint
        integer                     :: idomain_g, ielement_g, mapping, nnodes_element
        integer                     :: max_element_nodes, nnodes
        integer,    allocatable     :: nodes(:)
        integer(ik), allocatable    :: conn(:,:)



        ! Receive number of connectivities and inititalize partition
        call MPI_Recv(nconn,1,MPI_INTEGER, GLOBAL_MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        call partition%init(nconn)


        ! Receive connectivities
        do iconn = 1,nconn


            ! Recv connectivity information
            call MPI_Recv(dims,3,MPI_INTEGER, GLOBAL_MASTER, 1+iconn, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            nelements         = dims(1)
            max_element_nodes = dims(2)
            nnodes            = dims(3)

            call partition%connectivities(iconn)%init(nelements,nnodes)

            if (allocated(conn)) deallocate(conn)
            allocate(conn(nelements,max_element_nodes), stat=ierr)
            if (ierr /= 0) call AllocationError



            ! Recv connectivity data
            call MPI_Recv(conn,dims(1)*dims(2), MPI_INTEGER4, GLOBAL_MASTER, 1+nconn+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)




            ! Assemble connectivity
            do ielem = 1,nelements
                idomain_g  = conn(ielem,1)
                ielement_g = conn(ielem,2)
                mapping    = conn(ielem,3)
                nnodes_element = (mapping+1)*(mapping+1)*(mapping+1)
                nodes          = conn(ielem,4:4+nnodes_element-1)

                call partition%connectivities(iconn)%data(ielem)%init(mapping)
                call partition%connectivities(iconn)%data(ielem)%set_domain_index(idomain_g)
                call partition%connectivities(iconn)%data(ielem)%set_element_index(ielement_g)
                call partition%connectivities(iconn)%data(ielem)%set_element_mapping(mapping)
                call partition%connectivities(iconn)%data(ielem)%set_element_nodes(nodes)

            end do

        end do !iconn





        
!        do iprint = 0,NRANK-1
!            if (iprint == IRANK) then
!                print*, 'Hi from rank: ', IRANK
!                do ielem = 1,nelements
!                    print*, partition%connectivities(1)%data(ielem)%data
!                end do
!            end if
!            call MPI_Barrier(MPI_COMM_WORLD,ierr)
!        end do

    end subroutine recv_partition
    !***************************************************************************************************




















end module mod_partitioners
