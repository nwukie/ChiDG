module mod_partitioners
#include <messenger.h>
    use type_meshdata,  only: meshdata_t
    use mod_chidg_mpi,  only: NRANK, IRANK
    use mod_metis,      only: METIS_PartMeshNodal

    use type_connectivity,  only: connectivity_t
    use type_partition,     only: partition_t

    implicit none







contains




    !>  Take connectivity information from multiple domains, combine them into one
    !!  global connectivity, and partition.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/8/2016
    !! 
    !!
    !!  @param[in]      connectivities  Array of connectivity types, one for each domain
    !!  @param[inout]   partitions      Mesh partition data returned
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine partition_connectivity(connectivities, partitions)
        use iso_c_binding,  only: c_ptr, c_int, c_null_ptr
        type(connectivity_t),           intent(in)    :: connectivities(:)
        type(partition_t), allocatable, intent(inout) :: partitions(:)


        logical :: serial = .false.
        logical :: parallel = .false.
        logical :: partition_blocks = .false.


        integer                     :: idomain, ndomains, ielem_domain, nelem_domain, ierr
        integer                     :: ielem, inode, ielem_node, node_offset

        integer(c_int)              :: nelem, nnodes
        integer(c_int), allocatable :: eptr(:), eind(:), epart(:), npart(:)
        integer(c_int)              :: n, npartitions
        type(c_ptr)                 :: vwgt, vsize, options, tpwgts

        ! Unused METIS parameters
        vwgt    = c_null_ptr
        vsize   = c_null_ptr
        options = c_null_ptr
        tpwgts  = c_null_ptr


        !
        ! Get number of domains
        !
        ndomains = size(connectivities)


        !
        ! Accumulate total number of elements and nodes across all domains
        !
        nelem  = 0
        nnodes = 0
        do idomain = 1,ndomains
            nelem  = nelem  + connectivities(idomain)%get_nelements()
            nnodes = nnodes + connectivities(idomain)%get_nnodes()
        end do


        !
        ! Allocate storage for partitioning
        !
        allocate(eptr(nelem+1), eind(nelem*8), epart(nelem), npart(nnodes), partitions(NRANK), stat=ierr)
        if (ierr /= 0) call AllocationError





        !
        ! Assemble connectivities into METIS-readable format
        !
        inode = 1
        ielem = 1
        node_offset = 0
        do idomain = 1,ndomains

            nelem_domain = connectivities(idomain)%get_nelements()

            do ielem_domain = 1,nelem_domain

                ! Set pointer index to beginning of element connectivity
                eptr(ielem) = (ielem-1) * 8     ! 0-based


                ! Set element connectivity
                do ielem_node = 1,8
                    eind(inode) = connectivities(idomain)%data(ielem_domain,2+ielem_node) + node_offset - 1    ! 0-based
                    inode = inode + 1
                end do


                ielem = ielem + 1
            end do ! ielem_domain
        

            ! Offset by number of nodes in all previous domains
            node_offset = node_offset + connectivities(idomain)%get_nnodes()

        end do !idomain





        !
        ! METIS documentation is not clear why/what this last value is for.
        !
        eptr(nelem+1) = size(eind)



        !
        ! Partition mesh using METIS
        !
        npartitions = NRANK
        call METIS_PartMeshNodal(nelem,nnodes,eptr,eind,vwgt,vsize,npartitions,tpwgts,options,n,epart,npart)




        !
        ! Assemble partitions
        !
        do idomain = 1,ndomains





        end do ! idomain
        






    end subroutine partition_connectivity
    !*************************************************************************************************














    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/8/2016
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine send_partitions(partitions)
        type(partition_t),  intent(in)   :: partitions(:)


    end subroutine send_partitions
    !***************************************************************************************************



























end module mod_partitioners
