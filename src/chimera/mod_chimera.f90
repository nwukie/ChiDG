!> This module contains procedures for initializing and maintaining the Chimera
!! interfaces.
!!
!!  @author Nathan A. Wukie
!!  @date   2/1/2016
!!
!!
!---------------------------------------------------------------
module mod_chimera
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: NFACES, ORPHAN, CHIMERA, &
                                      X_DIR,  Y_DIR,   Z_DIR, &
                                      XI_DIR, ETA_DIR, ZETA_DIR, &
                                      ONE, ZERO, TWO, TWO_DIM, THREE_DIM, RKTOL, &
                                      INVALID_POINT, VALID_POINT, NO_PROC

    use type_mesh,              only: mesh_t
    use type_point,             only: point_t
    use type_element_info,      only: element_info_t
    use type_face_info,         only: face_info_t
    use type_ivector,           only: ivector_t
    use type_rvector,           only: rvector_t
    use type_pvector,           only: pvector_t

    use mod_polynomial,         only: polynomialVal
    use mod_periodic,           only: compute_periodic_offset
    use mod_chidg_mpi,          only: IRANK, NRANK, ChiDG_COMM
    use mpi_f08,                only: MPI_BCast, MPI_Send, MPI_Recv, MPI_INTEGER4, MPI_REAL8, MPI_LOGICAL, MPI_ANY_TAG, MPI_STATUS_IGNORE
    implicit none









contains


    !>  Routine for detecting Chimera faces. 
    !!
    !!  Routine flags face as a Chimera face if it has an ftype==ORPHAN, indicating it is not an interior
    !!  face and it has not been assigned a boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[inout]   mesh    Array of mesh types. One for each domain.
    !!
    !-----------------------------------------------------------------------------------------------------------------------
    subroutine detect_chimera_faces(mesh)
        type(mesh_t),   intent(inout)   :: mesh(:)

        integer(ik) :: idom, ndom, ielem, iface, ierr, nchimera_faces, ChiID
        logical     :: orphan_face = .false.
        logical     :: chimera_face = .false.

        !
        ! Get number of domains
        !
        ndom = size(mesh)

        
        !
        ! Loop through each element of each domain and look for ORPHAN face-types.
        ! If orphan is found, designate as CHIMERA and increment nchimera_faces
        !
        do idom = 1,ndom
            nchimera_faces = 0

            do ielem = 1,mesh(idom)%nelem


                !
                ! Loop through each face
                !
                do iface = 1,NFACES

                    !
                    ! Test if the current face is unattached
                    !
                    orphan_face = ( mesh(idom)%faces(ielem,iface)%ftype == ORPHAN ) 


                    !
                    ! If orphan_face, set as Chimera face so it can search for donors in other domains
                    !
                    if (orphan_face) then
                        ! Increment domain-local chimera face count
                        nchimera_faces = nchimera_faces + 1

                        ! Set face-type to CHIMERA
                        mesh(idom)%faces(ielem,iface)%ftype = CHIMERA

                        ! Set domain-local Chimera identifier. Really, just the index order which they were detected in, starting from 1.
                        ! The n-th chimera face
                        mesh(idom)%faces(ielem,iface)%ChiID = nchimera_faces
                    end if


                end do ! iface

            end do ! ielem



            !
            ! Set total number of Chimera faces detected for domain - idom
            !
            mesh(idom)%chimera%recv%nfaces = nchimera_faces


            !
            ! Allocate chimera_receiver_data for each chimera face in the current domain
            !
            allocate(mesh(idom)%chimera%recv%data(nchimera_faces), stat=ierr)
            if (ierr /= 0) call AllocationError


        end do ! idom






        !
        ! Now that all CHIMERA faces have been identified and we know the total number,
        ! we can store their data in the mesh-local chimera data container.
        !
        do idom = 1,ndom
            do ielem = 1,mesh(idom)%nelem

                !
                ! Loop through each face of current element
                !
                do iface = 1,NFACES

                    chimera_face = ( mesh(idom)%faces(ielem,iface)%ftype == CHIMERA )
                    if ( chimera_face ) then

                        !
                        ! Set receiver information for Chimera face
                        !
                        ChiID = mesh(idom)%faces(ielem,iface)%ChiID
                        mesh(idom)%chimera%recv%data(ChiID)%receiver_proc      = IRANK
                        mesh(idom)%chimera%recv%data(ChiID)%receiver_domain_g  = mesh(idom)%idomain_g
                        mesh(idom)%chimera%recv%data(ChiID)%receiver_domain_l  = mesh(idom)%idomain_l
                        mesh(idom)%chimera%recv%data(ChiID)%receiver_element_g = mesh(idom)%elems(ielem)%ielement_g
                        mesh(idom)%chimera%recv%data(ChiID)%receiver_element_l = mesh(idom)%elems(ielem)%ielement_l
                        mesh(idom)%chimera%recv%data(ChiID)%receiver_face      = iface
                    end if

                end do ! iface

            end do ! ielem
        end do ! idom



    end subroutine detect_chimera_faces
    !*********************************************************************************************************************















    !>  Routine for generating the data in a chimera_receiver_data instance. This includes donor_domain
    !!  and donor_element indices.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @parma[in]  mesh    Array of mesh_t instances
    !!
    !---------------------------------------------------------------------------------------------------------------------
    subroutine detect_chimera_donors(mesh)
        type(mesh_t),   intent(inout)   :: mesh(:)

        integer(ik) :: idom, igq, ichimera_face, idonor, ierr, iproc, idonor_proc, iproc_loop
        integer(ik) :: ndonors, neqns, nterms_s
        integer(ik) :: idonor_domain_g, idonor_element_g
        integer(ik) :: idonor_domain_l, idonor_element_l
        integer(ik) :: idomain_g_list, idomain_l_list, ielement_g_list, ielement_l_list, neqns_list, nterms_s_list, nterms_c_list, iproc_list
        integer(ik) :: local_domain_g, parallel_domain_g, donor_domain_g, donor_index
        integer(ik), allocatable    :: domains_g(:)

        integer(ik) :: receiver_indices(5), parallel_indices(8)
        real(rk)    :: gq_coords(3), parallel_coords(3), donor_vol, local_vol, parallel_vol
        real(rk), allocatable   :: donor_vols(:)

        real(rk)    :: offset_x, offset_y, offset_z

        type(face_info_t)           :: receiver
        type(element_info_t)        :: donor
        type(point_t)               :: donor_coord
        type(point_t)               :: gq_node
        type(point_t)               :: dummy_coord
        logical                     :: new_donor     = .false.
        logical                     :: already_added = .false.
        logical                     :: donor_match   = .false.
        logical                     :: searching
        logical                     :: donor_found
        logical                     :: proc_has_donor
        logical                     :: still_need_donor
        logical                     :: local_donor, parallel_donor
        logical                     :: use_local, use_parallel, get_donor

        type(ivector_t)             :: ddomain_g, ddomain_l, delement_g, delement_l, dproc, dneqns, dnterms_s, dnterms_c, donor_procs
        type(ivector_t)             :: donor_proc_indices, donor_proc_domains
        type(rvector_t)             :: donor_proc_vols
        type(pvector_t)             :: dcoordinate



        !
        ! Loop through processes. One will process its chimera faces and try to find processor-local donors. If it can't
        ! find on-processor donors, then it will broadcast a search request to all other processors. All other processors
        ! receive the request and return if they have a donor element or not.
        !
        ! The processes loop through in this serial fashion until all processors have processed their Chimera faces and
        ! have found donor elements for the quadrature nodes.
        !
        do iproc = 0,NRANK-1






            !
            ! iproc searches for donors for it's Chimera faces
            !
            if ( iproc == IRANK ) then
    

                !
                ! Loop over domains
                !
                do idom = 1,size(mesh)

                    call write_line('Detecting chimera donors for domain: ', idom, delimiter='  ')


                    !
                    ! Loop over faces and process Chimera-type faces
                    !
                    do ichimera_face = 1,mesh(idom)%chimera%recv%nfaces

                        !
                        ! Get location of the face receiving Chimera data
                        !
                        receiver%idomain_g  = mesh(idom)%chimera%recv%data(ichimera_face)%receiver_domain_g
                        receiver%idomain_l  = mesh(idom)%chimera%recv%data(ichimera_face)%receiver_domain_l
                        receiver%ielement_g = mesh(idom)%chimera%recv%data(ichimera_face)%receiver_element_g
                        receiver%ielement_l = mesh(idom)%chimera%recv%data(ichimera_face)%receiver_element_l
                        receiver%iface      = mesh(idom)%chimera%recv%data(ichimera_face)%receiver_face

                        call write_line('   Face ', ichimera_face,' of ',mesh(idom)%chimera%recv%nfaces, delimiter='  ')

                        !
                        ! Loop through quadrature nodes on Chimera face and find donors
                        !
                        do igq = 1,mesh(receiver%idomain_l)%faces(receiver%ielement_l,receiver%iface)%gq%face%nnodes

                            !
                            ! Get node coordinates
                            !
                            gq_node = mesh(receiver%idomain_l)%faces(receiver%ielement_l,receiver%iface)%quad_pts(igq)


                            !
                            ! Get offset coordinates from face for potential periodic offset.
                            !
                            call compute_periodic_offset(mesh(receiver%idomain_l)%faces(receiver%ielement_l,receiver%iface), gq_node, offset_x, offset_y, offset_z)

                            call gq_node%add_x(offset_x)
                            call gq_node%add_y(offset_y)
                            call gq_node%add_z(offset_z)


                            searching = .true.
                            call MPI_BCast(searching,1,MPI_LOGICAL, iproc, ChiDG_COMM, ierr)

                            ! Send gq node physical coordinates
                            gq_coords(1) = gq_node%c1_
                            gq_coords(2) = gq_node%c2_
                            gq_coords(3) = gq_node%c3_
                            call MPI_BCast(gq_coords,3,MPI_REAL8, iproc, ChiDG_COMM, ierr)

                            ! Send receiver indices
                            receiver_indices(1) = receiver%idomain_g
                            receiver_indices(2) = receiver%idomain_l
                            receiver_indices(3) = receiver%ielement_g
                            receiver_indices(4) = receiver%ielement_l
                            receiver_indices(5) = receiver%iface
                            call MPI_BCast(receiver_indices,5,MPI_INTEGER4, iproc, ChiDG_COMM, ierr)



                            !
                            ! Call routine to find LOCAL gq donor for current node
                            !
                            call find_gq_donor(mesh,gq_node, receiver, donor, donor_coord, donor_volume=local_vol)
                            donor_found = (donor_coord%status == VALID_POINT)

                            local_domain_g = 0
                            local_donor = .false.
                            if ( donor_found ) then
                                local_donor = .true.
                            end if



                            !
                            ! Check which processors have a valid donor
                            !
                            do idonor_proc = 0,NRANK-1
                                if (idonor_proc /= IRANK) then
                                    ! Check if a donor was found on proc idonor_proc
                                    call MPI_Recv(proc_has_donor,1,MPI_LOGICAL, idonor_proc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)

                                    if (proc_has_donor) then
                                        call donor_proc_indices%push_back(idonor_proc)

                                        call MPI_Recv(donor_domain_g,1,MPI_INTEGER4, idonor_proc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
                                        call donor_proc_domains%push_back(donor_domain_g)

                                        call MPI_Recv(donor_vol,1,MPI_REAL8, idonor_proc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
                                        call donor_proc_vols%push_back(donor_vol)
                                    end if
                                end if
                            end do !idonor_proc



                            !
                            !! If there is a parallel donor, determine which has the lowest domain index
                            ! If there is a parallel donor, determine which has the lowest volume
                            !
                            parallel_domain_g = 0
                            parallel_donor = .false.
                            if ( donor_proc_indices%size() > 0 ) then
                                donor_vols = donor_proc_vols%data()
                                donor_index = minloc(donor_vols,1)

                                parallel_vol = donor_vols(donor_index)
                                parallel_donor = .true.
                            end if
                            


                            !
                            ! Determine which donor to use
                            !
                            if ( local_donor .and. parallel_donor ) then
                                use_local    = (local_vol    < parallel_vol)
                                use_parallel = (parallel_vol < local_vol   )

                            else if (local_donor .and. (.not. parallel_donor)) then
                                use_local = .true.
                                use_parallel = .false.

                            else if (parallel_donor .and. (.not. local_donor)) then
                                use_local = .false.
                                use_parallel = .true.

                            else
                                call chidg_signal(FATAL,"detect_chimera_donor: no valid donor found")
                            end if




                            !
                            ! Overwrite donor data if use_parallel
                            !
                            if (use_parallel) then

                                !
                                ! Send a message to the processes that have donors to indicate if we would like to use it
                                !
                                do iproc_loop = 1,donor_proc_indices%size()
                                    idonor_proc = donor_proc_indices%at(iproc_loop)

                                    get_donor = (iproc_loop == donor_index)

                                    call MPI_Send(get_donor,1,MPI_LOGICAL, idonor_proc, 0, ChiDG_COMM, ierr)
                                end do !idonor_proc

                                !
                                ! Receive parallel donor index from processor indicated
                                !
                                idonor_proc = donor_proc_indices%at(donor_index)
                                call MPI_Recv(parallel_indices,8,MPI_INTEGER4, idonor_proc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
                                donor%idomain_g  = parallel_indices(1)
                                donor%idomain_l  = parallel_indices(2)
                                donor%ielement_g = parallel_indices(3)
                                donor%ielement_l = parallel_indices(4)
                                donor%iproc      = parallel_indices(5)
                                donor%neqns      = parallel_indices(6)
                                donor%nterms_s   = parallel_indices(7)
                                donor%nterms_c   = parallel_indices(8)



                                ! Receive donor local coordinate
                                call MPI_Recv(parallel_coords,3,MPI_REAL8, idonor_proc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
                                call donor_coord%set(parallel_coords(1), parallel_coords(2), parallel_coords(3))

                            else if (use_local) then

                                ! Send a message to all procs with donors that we don't need them
                                get_donor = .false.
                                do iproc_loop = 1,donor_proc_indices%size()
                                    idonor_proc = donor_proc_indices%at(iproc_loop)

                                    call MPI_Send(get_donor,1,MPI_LOGICAL, idonor_proc, 0, ChiDG_COMM, ierr)
                                end do

                            else 
                                call chidg_signal(FATAL,"detect_chimera_donors: no local or parallel donor found")

                            end if




                            !
                            ! Add donor location and coordinate
                            !
                            call ddomain_g%push_back(donor%idomain_g)
                            call ddomain_l%push_back(donor%idomain_l)
                            call delement_g%push_back(donor%ielement_g)
                            call delement_l%push_back(donor%ielement_l)
                            call dproc%push_back(donor%iproc)
                            call dneqns%push_back(donor%neqns)
                            call dnterms_s%push_back(donor%nterms_s)
                            call dnterms_c%push_back(donor%nterms_c)

                            call dcoordinate%push_back(donor_coord)


                            !
                            ! Clear working vectors
                            !
                            call donor_proc_indices%clear()
                            call donor_proc_domains%clear()
                            call donor_proc_vols%clear()

                        end do ! igq




                        !
                        ! Count number of unique donors to the current face and 
                        ! add to chimera donor data 
                        !
                        ndonors = 0
                        do igq = 1,ddomain_g%size()

                            idomain_g_list  = ddomain_g%at(igq)
                            idomain_l_list  = ddomain_l%at(igq)
                            ielement_g_list = delement_g%at(igq)
                            ielement_l_list = delement_l%at(igq)
                            iproc_list      = dproc%at(igq)
                            neqns_list      = dneqns%at(igq)
                            nterms_s_list   = dnterms_s%at(igq)
                            nterms_c_list   = dnterms_c%at(igq)


                            !
                            ! Check if domain/element pair has already been added to the chimera donor data
                            !
                            already_added = .false.
                            do idonor = 1,mesh(idom)%chimera%recv%data(ichimera_face)%donor_domain_g%size()

                                idonor_domain_g  = mesh(idom)%chimera%recv%data(ichimera_face)%donor_domain_g%at(idonor)
                                idonor_domain_l  = mesh(idom)%chimera%recv%data(ichimera_face)%donor_domain_l%at(idonor)
                                idonor_element_g = mesh(idom)%chimera%recv%data(ichimera_face)%donor_element_g%at(idonor)
                                idonor_element_l = mesh(idom)%chimera%recv%data(ichimera_face)%donor_element_l%at(idonor)
                                idonor_proc      = mesh(idom)%chimera%recv%data(ichimera_face)%donor_proc%at(idonor)

                                already_added = ( (idomain_g_list == idonor_domain_g)   .and. &
                                                  (idomain_l_list == idonor_domain_l)   .and. & 
                                                  (ielement_g_list == idonor_element_g) .and. &
                                                  (ielement_l_list == idonor_element_l)       &
                                                )
                                if (already_added) exit
                            end do
                            

                            !
                            ! If the current domain/element pair was not found in the chimera donor data, then add it
                            !
                            if (.not. already_added) then
                                call mesh(idom)%chimera%recv%data(ichimera_face)%donor_domain_g%push_back(idomain_g_list)
                                call mesh(idom)%chimera%recv%data(ichimera_face)%donor_domain_l%push_back(idomain_l_list)
                                call mesh(idom)%chimera%recv%data(ichimera_face)%donor_element_g%push_back(ielement_g_list)
                                call mesh(idom)%chimera%recv%data(ichimera_face)%donor_element_l%push_back(ielement_l_list)

                                call mesh(idom)%chimera%recv%data(ichimera_face)%donor_proc%push_back(iproc_list)
                                call mesh(idom)%chimera%recv%data(ichimera_face)%donor_neqns%push_back(neqns_list)
                                call mesh(idom)%chimera%recv%data(ichimera_face)%donor_nterms_s%push_back(nterms_s_list)

                                ! Initialize storage
                                call mesh(idom)%chimera%recv%data(ichimera_face)%donor_recv_comm%push_back(0)
                                call mesh(idom)%chimera%recv%data(ichimera_face)%donor_recv_domain%push_back(0)
                                call mesh(idom)%chimera%recv%data(ichimera_face)%donor_recv_element%push_back(0)
                                ndonors = ndonors + 1
                            end if

                        end do ! igq



                        !
                        ! Allocate chimera donor coordinate and quadrature index arrays. One list for each donor
                        !
                        allocate( mesh(idom)%chimera%recv%data(ichimera_face)%donor_coords(ndonors), &
                                  mesh(idom)%chimera%recv%data(ichimera_face)%donor_gq_indices(ndonors), stat=ierr)
                        if (ierr /= 0) call AllocationError





                        !
                        ! Now save donor coordinates and gq indices to their appropriate donor list
                        !
                        do igq = 1,ddomain_g%size()


                            idomain_g_list  = ddomain_g%at(igq)
                            idomain_l_list  = ddomain_l%at(igq)
                            ielement_g_list = delement_g%at(igq)
                            ielement_l_list = delement_l%at(igq)


                            !
                            ! Check if domain/element pair has already been added to the chimera donor data
                            !
                            donor_match = .false.
                            do idonor = 1,mesh(idom)%chimera%recv%data(ichimera_face)%donor_domain_g%size()

                                idonor_domain_g  = mesh(idom)%chimera%recv%data(ichimera_face)%donor_domain_g%at(idonor)
                                idonor_domain_l  = mesh(idom)%chimera%recv%data(ichimera_face)%donor_domain_l%at(idonor)
                                idonor_element_g = mesh(idom)%chimera%recv%data(ichimera_face)%donor_element_g%at(idonor)
                                idonor_element_l = mesh(idom)%chimera%recv%data(ichimera_face)%donor_element_l%at(idonor)

                                donor_match = ( (idomain_g_list == idonor_domain_g)   .and. &
                                                (idomain_l_list == idonor_domain_l)   .and. & 
                                                (ielement_g_list == idonor_element_g) .and. &
                                                (ielement_l_list == idonor_element_l) )

                                if (donor_match) then
                                    call mesh(idom)%chimera%recv%data(ichimera_face)%donor_gq_indices(idonor)%push_back(igq)
                                    call mesh(idom)%chimera%recv%data(ichimera_face)%donor_coords(idonor)%push_back(dcoordinate%at(igq))
                                    exit
                                end if
                            end do
                            
                        end do







                        !
                        ! Clear temporary donor arrays
                        !
                        call ddomain_g%clear()
                        call ddomain_l%clear()
                        call delement_g%clear()
                        call delement_l%clear()
                        call dcoordinate%clear()
                        call dproc%clear()
                        call dneqns%clear()
                        call dnterms_s%clear()
                        call dnterms_c%clear()



                    end do ! iface

                end do ! idom

            searching = .false.
            call MPI_BCast(searching,1,MPI_LOGICAL, iproc, ChiDG_COMM, ierr)

            end if ! iproc == IRANK










            !
            ! Each other proc waits for donor requests from iproc and sends back donors if they are found.
            !
            if (iproc /= IRANK) then

                !
                ! Check if iproc is searching for a node 
                !
                call MPI_BCast(searching,1,MPI_LOGICAL,iproc,ChiDG_COMM,ierr)



                do while(searching)

                    !
                    ! Receive gq node physical coordinates from iproc
                    !
                    call MPI_BCast(gq_coords,3,MPI_REAL8, iproc, ChiDG_COMM, ierr)
                    call gq_node%set(gq_coords(1), gq_coords(2), gq_coords(3))

                    
                    !
                    ! Receive receiver indices
                    !
                    call MPI_BCast(receiver_indices,5,MPI_INTEGER4, iproc, ChiDG_COMM, ierr)
                    receiver%idomain_g  = receiver_indices(1)
                    receiver%idomain_l  = receiver_indices(2)
                    receiver%ielement_g = receiver_indices(3)
                    receiver%ielement_l = receiver_indices(4)
                    receiver%iface      = receiver_indices(5)


                    !
                    ! Try to find donor
                    !
                    call find_gq_donor(mesh,gq_node,receiver,donor,donor_coord, donor_volume=donor_vol)
                    donor_found = (donor_coord%status == VALID_POINT)

                    
                    !
                    ! Send status
                    !
                    call MPI_Send(donor_found,1,MPI_LOGICAL,iproc,0,ChiDG_COMM,ierr)

                    if (donor_found) then

                        call MPI_Send(donor%idomain_g,1,MPI_INTEGER4,iproc,0,ChiDG_COMM,ierr)
                        call MPI_Send(donor_vol,1,MPI_REAL8,iproc,0,ChiDG_COMM,ierr)

                        call MPI_Recv(still_need_donor,1,MPI_LOGICAL, iproc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)

                        if (still_need_donor) then
                            ! Add donor to the send list
                            call mesh(donor%idomain_l)%chimera%send%add_donor(donor%idomain_g, donor%idomain_l, donor%ielement_g, donor%ielement_l, iproc)


                            ! Send donor indices
                            parallel_indices(1) = donor%idomain_g
                            parallel_indices(2) = donor%idomain_l
                            parallel_indices(3) = donor%ielement_g
                            parallel_indices(4) = donor%ielement_l
                            parallel_indices(5) = donor%iproc
                            parallel_indices(6) = donor%neqns
                            parallel_indices(7) = donor%nterms_s
                            parallel_indices(8) = donor%nterms_c
                            call MPI_Send(parallel_indices,8,MPI_INTEGER4,iproc,0,ChiDG_COMM,ierr)

                            ! Send donor-local coordinate for the quadrature node
                            parallel_coords(1) = donor_coord%c1_
                            parallel_coords(2) = donor_coord%c2_
                            parallel_coords(3) = donor_coord%c3_
                            call MPI_Send(parallel_coords,3,MPI_REAL8,iproc,0,ChiDG_COMM,ierr)
                        end if

                    end if


                    !
                    ! Check if iproc is searching for another node
                    !
                    call MPI_BCast(searching,1,MPI_LOGICAL,iproc,ChiDG_COMM,ierr)


                end do ! while searching



            end if ! iproc /= IRANK




        end do ! iproc




    end subroutine detect_chimera_donors
    !***********************************************************************************************************************














    !>  Find the domain and element indices for an element that contains a given quadrature node and can donate
    !!  interpolated solution values to the receiver face.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]      mesh            Array of mesh_t instances
    !!  @param[in]      gq_node         GQ point that needs to find a donor
    !!  @param[in]      receiver_face   Location of face containing the gq_node
    !!  @param[inout]   donor_element   Location of the donor element that was found
    !!  @param[inout]   donor_coord     Point defining the location of the GQ point in the donor coordinate system
    !!
    !-----------------------------------------------------------------------------------------------------------------------
    subroutine find_gq_donor(mesh,gq_node,receiver_face,donor_element,donor_coordinate,donor_volume)
        type(mesh_t),               intent(in)              :: mesh(:)
        type(point_t),              intent(in)              :: gq_node
        type(face_info_t),          intent(in)              :: receiver_face
        type(element_info_t),       intent(inout)           :: donor_element
        type(point_t),              intent(inout)           :: donor_coordinate
        real(rk),                   intent(inout), optional :: donor_volume


        integer(ik)             :: idom, ielem, inewton, spacedim
        integer(ik)             :: idomain_g, idomain_l, ielement_g, ielement_l
        integer(ik)             :: icandidate, ncandidates, idonor, ndonors
        integer(ik), allocatable    :: domains_g(:)
        type(point_t)           :: gq_comp
        real(rk)                :: xgq, ygq, zgq
        real(rk)                :: dx, dy, dz
        real(rk)                :: xi,  eta, zeta
        real(rk)                :: xn,  yn,  zn
        real(rk)                :: xmin, xmax, ymin, ymax, zmin, zmax
        real(rk)                :: tol
        type(ivector_t)         :: candidate_domains_g, candidate_domains_l, candidate_elements_g, candidate_elements_l
        type(ivector_t)         :: donors
        type(rvector_t)         :: donors_xi, donors_eta, donors_zeta

        integer(ik)             :: donor_index
        real(rk), allocatable   :: xcenter(:), ycenter(:), zcenter(:), dist(:), donor_vols(:)
        real(rk)                :: xcenter_recv, ycenter_recv, zcenter_recv

        logical                 :: contained = .false.
        logical                 :: receiver  = .false.
        logical                 :: donor_found



        tol = 10._rk*RKTOL


        xgq = gq_node%c1_
        ygq = gq_node%c2_
        zgq = gq_node%c3_


        !
        ! Loop through LOCAL domains and search for potential donor candidates
        !
        ncandidates = 0
        do idom = 1,size(mesh)
            idomain_g = mesh(idom)%idomain_g
            idomain_l = mesh(idom)%idomain_l



            !
            ! Loop through elements in the current domain
            !
            do ielem = 1,mesh(idom)%nelem
                ielement_g = mesh(idom)%elems(ielem)%ielement_g
                ielement_l = mesh(idom)%elems(ielem)%ielement_l

                !
                ! Get bounding coordinates for the current element
                !
                xmin = minval(mesh(idom)%elems(ielem)%elem_pts(:)%c1_)
                xmax = maxval(mesh(idom)%elems(ielem)%elem_pts(:)%c1_)
                ymin = minval(mesh(idom)%elems(ielem)%elem_pts(:)%c2_)
                ymax = maxval(mesh(idom)%elems(ielem)%elem_pts(:)%c2_)
                zmin = minval(mesh(idom)%elems(ielem)%elem_pts(:)%c3_)
                zmax = maxval(mesh(idom)%elems(ielem)%elem_pts(:)%c3_)

                !
                ! Grow bounding box by 10%. Use delta x,y,z instead of scaling xmin etc. in case xmin is 0
                !
                dx = abs(xmax - xmin)  
                dy = abs(ymax - ymin)
                dz = abs(zmax - zmin)

                xmin = xmin - 0.1*dx
                xmax = xmax + 0.1*dx
                ymin = ymin - 0.1*dy
                ymax = ymax + 0.1*dy
                zmin = (zmin-0.001) - 0.1*dz    ! This is to help 2D
                zmax = (zmax+0.001) + 0.1*dz    ! This is to help 2D

                !
                ! Test if gq_node is contained within the bounding coordinates
                !
                contained = ( (xmin < xgq) .and. (xgq < xmax ) .and. &
                              (ymin < ygq) .and. (ygq < ymax ) .and. &
                              (zmin < zgq) .and. (zgq < zmax ) )

                !
                ! Make sure that we arent adding the receiver element itself as a potential donor
                !
                receiver = ( (idomain_g == receiver_face%idomain_g) .and. (ielement_g == receiver_face%ielement_g) )


                !
                ! If the node was within the bounding coordinates, flag the element as a potential donor
                !
                if (contained .and. (.not. receiver)) then
                    call candidate_domains_g%push_back(idomain_g)
                    call candidate_domains_l%push_back(idomain_l)
                    call candidate_elements_g%push_back(ielement_g)
                    call candidate_elements_l%push_back(ielement_l)
                    ncandidates = ncandidates + 1
                end if


            end do ! ielem

        end do ! idom





        !
        ! Test gq_node physical coordinates on candidate element volume to try and map to donor local coordinates
        !
        ndonors = 0
        do icandidate = 1,ncandidates
            idomain_g  = candidate_domains_g%at(icandidate)
            idomain_l  = candidate_domains_l%at(icandidate)
            ielement_g = candidate_elements_g%at(icandidate)
            ielement_l = candidate_elements_l%at(icandidate)

            
            !
            ! Try to find donor (xi,eta,zeta) coordinates for receiver (xgq,ygq,zgq)
            !
            gq_comp = mesh(idomain_l)%elems(ielement_l)%computational_point(xgq,ygq,zgq)    ! Newton's method routine



            ! Add donor if gq_comp point is valid
            !
            ! NOTE: the exit call here enforces that only one candidate is considered and others are thrown away. Maybe not the best choice.
            !
            donor_found = (gq_comp%status == VALID_POINT)
            if ( donor_found ) then
                ndonors = ndonors + 1
                call donors%push_back(icandidate)
                call donors_xi%push_back(  gq_comp%c1_)
                call donors_eta%push_back( gq_comp%c2_)
                call donors_zeta%push_back(gq_comp%c3_)
                !exit
            end if

        end do ! icandidate






        




        !
        ! Sanity check on donors and set donor_element location
        !
        if (ndonors == 0) then
            donor_element%idomain_g  = 0
            donor_element%idomain_l  = 0
            donor_element%ielement_g = 0
            donor_element%ielement_l = 0
            donor_element%iproc      = NO_PROC

            donor_coordinate%status = INVALID_POINT



        elseif (ndonors == 1) then
            idonor = donors%at(1)   ! donor index from candidates
            donor_element%idomain_g  = candidate_domains_g%at(idonor)
            donor_element%idomain_l  = candidate_domains_l%at(idonor)
            donor_element%ielement_g = candidate_elements_g%at(idonor)
            donor_element%ielement_l = candidate_elements_l%at(idonor)
            donor_element%iproc      = IRANK
            donor_element%neqns      = mesh(donor_element%idomain_l)%elems(donor_element%ielement_l)%neqns
            donor_element%nterms_s   = mesh(donor_element%idomain_l)%elems(donor_element%ielement_l)%nterms_s
            donor_element%nterms_c   = mesh(donor_element%idomain_l)%elems(donor_element%ielement_l)%nterms_c

            xi   = donors_xi%at(1)
            eta  = donors_eta%at(1)
            zeta = donors_zeta%at(1)
            call donor_coordinate%set(xi,eta,zeta)
            donor_coordinate%status = VALID_POINT
            if (present(donor_volume)) donor_volume = mesh(donor_element%idomain_l)%elems(donor_element%ielement_l)%vol



        elseif (ndonors > 1) then
            !TODO: Account for case of multiple overlapping donors. When a gq node could be filled by two or more elements.
            !      Maybe, just choose one. Maybe, average contribution from all potential donors.
            !call chidg_signal(FATAL,"find_gq_donor: Multiple donors found for the same gq_node")

            !!
            !! OPTION 1: Choose donor with minimum global domain index
            !!
!
!            ! Get domain global indices of valid donors
!            if (allocated(domains_g)) deallocate(domains_g)
!            allocate(domains_g(donors%size()))
!            do idonor = 1,donors%size()
!                domains_g(idonor) = candidate_domains_g%at(donors%at(idonor))
!            end do
!
!            ! Get index of minimum domain
!            donor_index = minloc(domains_g,1)
!            idonor = donors%at(donor_index)
!
!            donor_element%idomain_g  = candidate_domains_g%at(idonor)
!            donor_element%idomain_l  = candidate_domains_l%at(idonor)
!            donor_element%ielement_g = candidate_elements_g%at(idonor)
!            donor_element%ielement_l = candidate_elements_l%at(idonor)
!            donor_element%iproc      = IRANK
!            donor_element%neqns      = mesh(donor_element%idomain_l)%elems(donor_element%ielement_l)%neqns
!            donor_element%nterms_s   = mesh(donor_element%idomain_l)%elems(donor_element%ielement_l)%nterms_s
!            donor_element%nterms_c   = mesh(donor_element%idomain_l)%elems(donor_element%ielement_l)%nterms_c
!
!
!            !!
!            !! OPTION 2: Choose donor with cell center closer to receiver center
!            !!
!
!            ! Get donor cell centers
!            if (allocated(xcenter) .or. allocated(ycenter) .or. allocated(zcenter) .or. allocated(dist)) deallocate(xcenter, ycenter, zcenter, dist)
!            allocate(xcenter(donors%size()), ycenter(donors%size()), zcenter(donors%size()), dist(donors%size()))
!            
!            do idonor = 1,donors%size()
!                xcenter(idonor) = mesh(candidate_domains_l%at(donors%at(idonor)))%elems(candidate_elements_l%at(donors%at(idonor)))%coords%getterm(1,1)
!                ycenter(idonor) = mesh(candidate_domains_l%at(donors%at(idonor)))%elems(candidate_elements_l%at(donors%at(idonor)))%coords%getterm(2,1)
!                zcenter(idonor) = mesh(candidate_domains_l%at(donors%at(idonor)))%elems(candidate_elements_l%at(donors%at(idonor)))%coords%getterm(3,1)
!            end do 
!
!            ! Get receiver cell center
!            xcenter_recv = mesh(receiver_face%idomain_l)%elems(receiver_face%ielement_l)%coords%getterm(1,1)
!            ycenter_recv = mesh(receiver_face%idomain_l)%elems(receiver_face%ielement_l)%coords%getterm(2,1)
!            zcenter_recv = mesh(receiver_face%idomain_l)%elems(receiver_face%ielement_l)%coords%getterm(3,1)
!
!            ! Compute the distance between receiver and donor cell centers
!            do idonor = 1,donors%size()
!                dist(idonor) = sqrt((xcenter(idonor)-xcenter_recv)**TWO  +  (ycenter(idonor)-ycenter_recv)**TWO  +  (zcenter(idonor)-zcenter_recv)**TWO)
!            end do
!
!            ! Select the donor element with the minimum distance between cell centers
!            donor_index = minloc(dist,1)
!            idonor = donors%at(donor_index)
!
!            ! Store donor data to be returned
!            donor_element%idomain_g  = candidate_domains_g%at(idonor)
!            donor_element%idomain_l  = candidate_domains_l%at(idonor)
!            donor_element%ielement_g = candidate_elements_g%at(idonor)
!            donor_element%ielement_l = candidate_elements_l%at(idonor)
!            donor_element%iproc      = IRANK
!            donor_element%neqns      = mesh(donor_element%idomain_l)%elems(donor_element%ielement_l)%neqns
!            donor_element%nterms_s   = mesh(donor_element%idomain_l)%elems(donor_element%ielement_l)%nterms_s
!            donor_element%nterms_c   = mesh(donor_element%idomain_l)%elems(donor_element%ielement_l)%nterms_c
!

            !!
            !! OPTION 3: Choose donor with minimum volume: should be best resolved
            !!
            if (allocated(donor_vols) ) deallocate(donor_vols)
            allocate(donor_vols(donors%size()))
            
            do idonor = 1,donors%size()
                donor_vols(idonor) = mesh(candidate_domains_l%at(donors%at(idonor)))%elems(candidate_elements_l%at(donors%at(idonor)))%vol
            end do 
    

            ! Get index of minimum domain
            donor_index = minloc(donor_vols,1)
            idonor = donors%at(donor_index)

            donor_element%idomain_g  = candidate_domains_g%at(idonor)
            donor_element%idomain_l  = candidate_domains_l%at(idonor)
            donor_element%ielement_g = candidate_elements_g%at(idonor)
            donor_element%ielement_l = candidate_elements_l%at(idonor)
            donor_element%iproc      = IRANK
            donor_element%neqns      = mesh(donor_element%idomain_l)%elems(donor_element%ielement_l)%neqns
            donor_element%nterms_s   = mesh(donor_element%idomain_l)%elems(donor_element%ielement_l)%nterms_s
            donor_element%nterms_c   = mesh(donor_element%idomain_l)%elems(donor_element%ielement_l)%nterms_c





            !
            ! Set donor coordinate and volume if present
            !
            xi   = donors_xi%at(donor_index)
            eta  = donors_eta%at(donor_index)
            zeta = donors_zeta%at(donor_index)
            call donor_coordinate%set(xi,eta,zeta)
            donor_coordinate%status = VALID_POINT
            if (present(donor_volume)) donor_volume = mesh(donor_element%idomain_l)%elems(donor_element%ielement_l)%vol


        else
            call chidg_signal(FATAL,"find_gq_donor: invalid number of donors")
        end if




    end subroutine find_gq_donor
    !****************************************************************************************************************















    !> Compute the matrices that interpolate solution data from a donor element expansion
    !! to the receiver nodes.
    !!
    !! These matrices get stored in mesh(idom)%chimera%recv%data(ChiID)%donor_interpolator
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-----------------------------------------------------------------------------------------------------------------
    subroutine compute_chimera_interpolators(mesh)
        type(mesh_t),   intent(inout)   :: mesh(:)

        integer(ik) :: idom, iChiID, idonor, ierr, ipt, iterm
        integer(ik) :: donor_idomain_g, donor_idomain_l, donor_ielement_g, donor_ielement_l
        integer(ik) :: npts, donor_nterms_s, spacedim

        type(point_t)           :: node
        real(rk), allocatable   :: interpolator(:,:)



        !
        ! Loop over all domains
        !
        do idom = 1,size(mesh)

            spacedim = mesh(idom)%spacedim

            !
            ! Loop over each chimera face
            !
            do iChiID = 1,mesh(idom)%chimera%recv%nfaces

                
                !
                ! For each donor, compute an interpolation matrix
                !
                do idonor = 1,mesh(idom)%chimera%recv%data(iChiID)%ndonors()

                    donor_idomain_g  = mesh(idom)%chimera%recv%data(iChiID)%donor_domain_g%at(idonor)
                    donor_idomain_l  = mesh(idom)%chimera%recv%data(iChiID)%donor_domain_l%at(idonor)
                    donor_ielement_g = mesh(idom)%chimera%recv%data(iChiID)%donor_element_g%at(idonor)
                    donor_ielement_l = mesh(idom)%chimera%recv%data(iChiID)%donor_element_l%at(idonor)
                    donor_nterms_s   = mesh(idom)%chimera%recv%data(iChiID)%donor_nterms_s%at(idonor)

                    !
                    ! Get number of GQ points this donor is responsible for
                    !
                    npts   = mesh(idom)%chimera%recv%data(iChiID)%donor_coords(idonor)%size()

                    !
                    ! Allocate interpolator matrix
                    !
                    if (allocated(interpolator)) deallocate(interpolator)
                    allocate(interpolator(npts,donor_nterms_s), stat=ierr)
                    if (ierr /= 0) call AllocationError

                    !
                    ! Compute values of modal polynomials at the donor nodes
                    !
                    do iterm = 1,donor_nterms_s
                        do ipt = 1,npts

                            node = mesh(idom)%chimera%recv%data(iChiID)%donor_coords(idonor)%at(ipt)
                            interpolator(ipt,iterm) = polynomialVal(spacedim,donor_nterms_s,iterm,node)

                        end do ! ipt
                    end do ! iterm

                    !
                    ! Store interpolator
                    !
                    call mesh(idom)%chimera%recv%data(iChiID)%donor_interpolator%push_back(interpolator)

                end do  ! idonor



            end do  ! iChiID
        end do  ! idom


    end subroutine compute_chimera_interpolators
    !******************************************************************************************************************








end module mod_chimera
