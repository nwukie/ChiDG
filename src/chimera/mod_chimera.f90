!>  This module contains procedures for initializing and maintaining the Chimera
!!  interfaces.
!!
!!  detect_chimera_faces
!!  detect_chimera_donors
!!  compute_chimera_interpolators
!!  find_gq_donor
!!
!!  detect_chimera_faces, detect_chimera_donors, and compute_chimera_interpolators are probably
!!  called in src/parallel/mod_communication.f90%establish_chimera_communication.
!!
!!  @author Nathan A. Wukie
!!  @date   2/1/2016
!!
!!
!-----------------------------------------------------------------------------------------------------
module mod_chimera
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: NFACES, ORPHAN, CHIMERA, &
                                      XI_DIR, ETA_DIR, ZETA_DIR, &
                                      ONE, ZERO, TWO, TWO_DIM, THREE_DIM, &
                                      INVALID_POINT, VALID_POINT, NO_PROC, NO_ID
    use mod_reference_elements, only: ref_elems

    use type_point
    use type_mesh,              only: mesh_t
    use type_chimera_send,      only: chimera_send
    use type_element_info,      only: element_info_t
    use type_face_info,         only: face_info_t, face_info
    use type_ivector,           only: ivector_t
    use type_rvector,           only: rvector_t
    use type_pvector,           only: pvector_t
    use type_mvector,           only: mvector_t

    use mod_polynomial,         only: polynomial_val, dpolynomial_val
    use mod_periodic,           only: get_periodic_offset
    use mod_chidg_mpi,          only: IRANK, NRANK, ChiDG_COMM
    use mpi_f08,                only: MPI_BCast, MPI_Send, MPI_Recv, MPI_INTEGER4, MPI_REAL8, &
                                      MPI_LOGICAL, MPI_ANY_TAG, MPI_STATUS_IGNORE
    use ieee_arithmetic,        only: ieee_is_nan
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
        type(mesh_t),   intent(inout)   :: mesh

        integer(ik) :: idom, ndom, ielem, iface, nnodes, ierr, ChiID
        logical     :: orphan_face = .false.
        logical     :: chimera_face = .false.

        !
        ! Get number of domains
        !
        ndom = mesh%ndomains()

        
        !
        ! Loop through each element of each domain and look for ORPHAN face-types.
        ! If orphan is found, designate as CHIMERA 
        !
        do idom = 1,ndom
            do ielem = 1,mesh%domain(idom)%nelem
                do iface = 1,NFACES

                    !
                    ! Test if the current face is unattached. 
                    ! Test also if the current face is CHIMERA in case this is being 
                    ! called as a reinitialization procedure.
                    !
                    orphan_face = ( mesh%domain(idom)%faces(ielem,iface)%ftype == ORPHAN .or. &
                                    mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA )


                    !
                    ! If orphan_face, set as Chimera face so it can search for donors in other domains
                    !
                    if (orphan_face) then

                        ! Set face-type to CHIMERA
                        mesh%domain(idom)%faces(ielem,iface)%ftype = CHIMERA
                        nnodes = size(mesh%domain(idom)%faces(ielem,iface)%jinv)

                        ! Set domain-local Chimera identifier. Really, just the index order which they were detected in, starting from 1.
                        ! The n-th chimera face
                        mesh%domain(idom)%faces(ielem,iface)%ChiID = mesh%domain(idom)%chimera%add_receiver(mesh%domain(idom)%idomain_g,                &
                                                                                                            mesh%domain(idom)%idomain_l,                &
                                                                                                            mesh%domain(idom)%elems(ielem)%ielement_g,  &
                                                                                                            mesh%domain(idom)%elems(ielem)%ielement_l,  &
                                                                                                            iface,                                      &
                                                                                                            nnodes)
                    end if


                end do ! iface
            end do ! ielem
        end do ! idom


    end subroutine detect_chimera_faces
    !*********************************************************************************************************************















    !>  Routine for generating the data in a chimera_receiver_data instance. This includes donor_domain
    !!  and donor_element indices.
    !!
    !!  For each Chimera face, find a donor for each quadrature node on the face, for a given node, initialize information
    !!  about its donor.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @parma[in]  mesh    Array of mesh_t instances
    !!
    !---------------------------------------------------------------------------------------------------------------------
    subroutine detect_chimera_donors(mesh)
        type(mesh_t),   intent(inout)   :: mesh

        integer(ik) :: idom, igq, ichimera_face, idonor, ierr, iproc,                           &
                       idonor_proc, iproc_loop,                                                 &
                       ndonors, neqns, nterms_s,                                                &
                       idonor_domain_g, idonor_element_g, idonor_domain_l, idonor_element_l,    &
                       idomain_g_list, idomain_l_list, ielement_g_list, ielement_l_list,        &
                       neqns_list, nterms_s_list, nterms_c_list, iproc_list, eqn_ID_list,       &
                       local_domain_g, parallel_domain_g, donor_domain_g, donor_index, donor_ID, send_ID
        integer(ik), allocatable    :: domains_g(:)
        integer(ik)                 :: receiver_indices(5), parallel_indices(9)


        real(rk)                :: donor_metric(3,3), parallel_metric(3,3)
        real(rk), allocatable   :: donor_vols(:)
        real(rk)                :: gq_coords(3), offset(3), gq_node(3), &
                                   donor_jinv, donor_vol, local_vol, parallel_vol, parallel_jinv
        real(rk)                :: d1dxi, d1deta, d1dzeta, &
                                   d2dxi, d2deta, d2dzeta, &
                                   d3dxi, d3deta, d3dzeta

        type(mvector_t) :: dmetric
        type(rvector_t) :: djinv

        type(face_info_t)           :: receiver
        type(element_info_t)        :: donor
        real(rk)                    :: donor_coord(3)
        logical                     :: new_donor     = .false.
        logical                     :: already_added = .false.
        logical                     :: donor_match   = .false.
        logical                     :: searching
        logical                     :: donor_found
        logical                     :: proc_has_donor
        logical                     :: still_need_donor
        logical                     :: local_donor, parallel_donor
        logical                     :: use_local, use_parallel, get_donor

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
    

                do idom = 1,mesh%ndomains()

                    call write_line('   Detecting chimera donors for domain: ', idom, delimiter='  ', ltrim=.false.)


                    !
                    ! Loop over faces and process Chimera-type faces
                    !
                    do ichimera_face = 1,mesh%domain(idom)%chimera%nreceivers()

                        !
                        ! Get location of the face receiving Chimera data
                        !
                        receiver%idomain_g  = mesh%domain(idom)%chimera%recv(ichimera_face)%idomain_g
                        receiver%idomain_l  = mesh%domain(idom)%chimera%recv(ichimera_face)%idomain_l
                        receiver%ielement_g = mesh%domain(idom)%chimera%recv(ichimera_face)%ielement_g
                        receiver%ielement_l = mesh%domain(idom)%chimera%recv(ichimera_face)%ielement_l
                        receiver%iface      = mesh%domain(idom)%chimera%recv(ichimera_face)%iface

                        call write_line('   Face ', ichimera_face,' of ',mesh%domain(idom)%chimera%nreceivers(), delimiter='  ')

                        !
                        ! Loop through quadrature nodes on Chimera face and find donors
                        !
                        do igq = 1,mesh%domain(receiver%idomain_l)%faces(receiver%ielement_l,receiver%iface)%basis_s%nnodes_if()

                            !
                            ! Get node coordinates
                            !
                            gq_node = mesh%domain(receiver%idomain_l)%faces(receiver%ielement_l,receiver%iface)%quad_pts(igq,1:3)


                            !
                            ! Get offset coordinates from face for potential periodic offset.
                            !
                            offset = get_periodic_offset(mesh%domain(receiver%idomain_l)%faces(receiver%ielement_l,receiver%iface))

                            gq_node = gq_node + offset


                            searching = .true.
                            call MPI_BCast(searching,1,MPI_LOGICAL, iproc, ChiDG_COMM, ierr)

                            ! Send gq node physical coordinates
                            call MPI_BCast(gq_node,3,MPI_REAL8, iproc, ChiDG_COMM, ierr)

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
                            call find_gq_donor(mesh,point_t(gq_node), receiver, donor, donor_coord, donor_found, donor_volume=local_vol)

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
                                    get_donor   = (iproc_loop == donor_index)
                                    call MPI_Send(get_donor,1,MPI_LOGICAL, idonor_proc, 0, ChiDG_COMM, ierr)

                                end do !idonor_proc

                                !
                                ! Receive parallel donor index from processor indicated
                                !
                                idonor_proc = donor_proc_indices%at(donor_index)
                                call MPI_Recv(parallel_indices,9,MPI_INTEGER4, idonor_proc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
                                donor%idomain_g  = parallel_indices(1)
                                donor%idomain_l  = parallel_indices(2)
                                donor%ielement_g = parallel_indices(3)
                                donor%ielement_l = parallel_indices(4)
                                donor%iproc      = parallel_indices(5)
                                donor%eqn_ID     = parallel_indices(6)
                                donor%neqns      = parallel_indices(7)
                                donor%nterms_s   = parallel_indices(8)
                                donor%nterms_c   = parallel_indices(9)



                                ! 1: Receive donor local coordinate
                                ! 2: Receive donor metric matrix
                                ! 3: Receive donor inverse jacobian mapping
                                call MPI_Recv(donor_coord,3,MPI_REAL8, idonor_proc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
                                call MPI_Recv(donor_metric,9,MPI_REAL8, idonor_proc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)
                                call MPI_Recv(donor_jinv,1,MPI_REAL8, idonor_proc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)



                            else if (use_local) then

                                ! Send a message to all procs with donors that we don't need them
                                get_donor = .false.
                                do iproc_loop = 1,donor_proc_indices%size()
                                    idonor_proc = donor_proc_indices%at(iproc_loop)
                                    call MPI_Send(get_donor,1,MPI_LOGICAL, idonor_proc, 0, ChiDG_COMM, ierr)
                                end do



                                ! Compute local metric
                                d1dxi   = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(1,1,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)
                                d2dxi   = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(2,1,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)
                                d3dxi   = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(3,1,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)
                                d1deta  = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(1,2,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)
                                d2deta  = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(2,2,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)
                                d3deta  = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(3,2,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)
                                d1dzeta = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(1,3,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)
                                d2dzeta = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(2,3,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)
                                d3dzeta = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(3,3,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)

                                donor_jinv = (d1dxi*d2deta*d3dzeta - d1deta*d2dxi*d3dzeta - &
                                              d1dxi*d2dzeta*d3deta + d1dzeta*d2dxi*d3deta + &
                                              d1deta*d2dzeta*d3dxi - d1dzeta*d2deta*d3dxi)


                                donor_metric(1,1) = (d2deta*d3dzeta - d2dzeta*d3deta)
                                donor_metric(2,1) = (d2dzeta*d3dxi  - d2dxi*d3dzeta)
                                donor_metric(3,1) = (d2dxi*d3deta   - d2deta*d3dxi)
                                donor_metric(1,2) = (d1dzeta*d3deta - d1deta*d3dzeta)
                                donor_metric(2,2) = (d1dxi*d3dzeta  - d1dzeta*d3dxi)
                                donor_metric(3,2) = (d1deta*d3dxi   - d1dxi*d3deta)
                                donor_metric(1,3) = (d1deta*d2dzeta - d1dzeta*d2deta)
                                donor_metric(2,3) = (d1dzeta*d2dxi  - d1dxi*d2dzeta)
                                donor_metric(3,3) = (d1dxi*d2deta   - d1deta*d2dxi)

                                ! Complete definition of metric term by scaling by J
                                donor_metric = donor_metric/donor_jinv


                            else 
                                call chidg_signal(FATAL,"detect_chimera_donors: no local or parallel donor found")

                            end if

                            !
                            ! Add donor
                            !
                            donor_ID = mesh%domain(idom)%chimera%recv(ichimera_face)%add_donor(donor%idomain_g, donor%idomain_l, donor%ielement_g, donor%ielement_l, donor%iproc)
                            call mesh%domain(idom)%chimera%recv(ichimera_face)%donor(donor_ID)%set_properties(donor%nterms_c,donor%nterms_s,donor%neqns,donor%eqn_ID)
                            call mesh%domain(idom)%chimera%recv(ichimera_face)%donor(donor_ID)%add_node(igq,donor_coord,donor_metric,donor_jinv)


                            !
                            ! Clear working vectors
                            !
                            call donor_proc_indices%clear()
                            call donor_proc_domains%clear()
                            call donor_proc_vols%clear()

                        end do ! igq

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
                    call MPI_BCast(gq_node,3,MPI_REAL8, iproc, ChiDG_COMM, ierr)

                    
                    !
                    ! Receive receiver indices
                    !
                    call MPI_BCast(receiver_indices,5,MPI_INTEGER4, iproc, ChiDG_COMM, ierr)
                    receiver = face_info(receiver_indices(1),   &
                                         receiver_indices(2),   &
                                         receiver_indices(3),   &
                                         receiver_indices(4),   &
                                         receiver_indices(5)    &
                                         )


                    !
                    ! Try to find donor
                    !
                    call find_gq_donor(mesh,point_t(gq_node),receiver,donor,donor_coord, donor_found, donor_volume=donor_vol)

                    
                    !
                    ! Send status
                    !
                    call MPI_Send(donor_found,1,MPI_LOGICAL,iproc,0,ChiDG_COMM,ierr)

                    if (donor_found) then

                        call MPI_Send(donor%idomain_g,1,MPI_INTEGER4,iproc,0,ChiDG_COMM,ierr)
                        call MPI_Send(donor_vol,1,MPI_REAL8,iproc,0,ChiDG_COMM,ierr)

                        call MPI_Recv(still_need_donor,1,MPI_LOGICAL, iproc, MPI_ANY_TAG, ChiDG_COMM, MPI_STATUS_IGNORE, ierr)

                        if (still_need_donor) then

                            ! Add donor to the chimera send collection
                            send_ID = mesh%domain(donor%idomain_l)%chimera%find_send(donor%idomain_g, donor%ielement_g)
                            if (send_ID == NO_ID) send_ID = mesh%domain(donor%idomain_l)%chimera%new_send()
                            mesh%domain(donor%idomain_l)%chimera%send(send_ID) = chimera_send(donor%idomain_g, donor%idomain_l, donor%ielement_g, donor%ielement_l)
                            call mesh%domain(donor%idomain_l)%chimera%send(send_ID)%send_procs%push_back(iproc)



                            ! 1: Send donor indices
                            ! 2: Send donor-local coordinate for the quadrature node
                            parallel_indices(1) = donor%idomain_g
                            parallel_indices(2) = donor%idomain_l
                            parallel_indices(3) = donor%ielement_g
                            parallel_indices(4) = donor%ielement_l
                            parallel_indices(5) = donor%iproc
                            parallel_indices(6) = donor%eqn_ID
                            parallel_indices(7) = donor%neqns
                            parallel_indices(8) = donor%nterms_s
                            parallel_indices(9) = donor%nterms_c

                            call MPI_Send(parallel_indices,9,MPI_INTEGER4,iproc,0,ChiDG_COMM,ierr)
                            call MPI_Send(donor_coord,3,MPI_REAL8,iproc,0,ChiDG_COMM,ierr)


                            ! Compute metric terms for the point in the donor element
                            d1dxi   = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(1,1,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)
                            d2dxi   = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(2,1,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)
                            d3dxi   = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(3,1,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)
                            d1deta  = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(1,2,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)
                            d2deta  = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(2,2,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)
                            d3deta  = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(3,2,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)
                            d1dzeta = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(1,3,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)
                            d2dzeta = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(2,3,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)
                            d3dzeta = mesh%domain(donor%idomain_l)%elems(donor%ielement_l)%metric_point(3,3,donor_coord(1),donor_coord(2),donor_coord(3),scale=.true.)


                            ! Compute inverse element jacobian
                            parallel_jinv = (d1dxi*d2deta*d3dzeta - d1deta*d2dxi*d3dzeta - &
                                             d1dxi*d2dzeta*d3deta + d1dzeta*d2dxi*d3deta + &
                                             d1deta*d2dzeta*d3dxi - d1dzeta*d2deta*d3dxi)


                            parallel_metric(1,1) = (d2deta*d3dzeta - d2dzeta*d3deta)
                            parallel_metric(2,1) = (d2dzeta*d3dxi  - d2dxi*d3dzeta)
                            parallel_metric(3,1) = (d2dxi*d3deta   - d2deta*d3dxi)
                            parallel_metric(1,2) = (d1dzeta*d3deta - d1deta*d3dzeta)
                            parallel_metric(2,2) = (d1dxi*d3dzeta  - d1dzeta*d3dxi)
                            parallel_metric(3,2) = (d1deta*d3dxi   - d1dxi*d3deta)
                            parallel_metric(1,3) = (d1deta*d2dzeta - d1dzeta*d2deta)
                            parallel_metric(2,3) = (d1dzeta*d2dxi  - d1dxi*d2dzeta)
                            parallel_metric(3,3) = (d1dxi*d2deta   - d1deta*d2dxi)

                            ! Complete definition of metric by scaling by J
                            parallel_metric = parallel_metric/parallel_jinv

                            ! Communicate metric and jacobian 
                            call MPI_Send(parallel_metric,9,MPI_REAL8,iproc,0,ChiDG_COMM,ierr)
                            call MPI_Send(parallel_jinv,1,MPI_REAL8,iproc,0,ChiDG_COMM,ierr)

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











    !> Compute the matrices that interpolate solution data from a donor element expansion
    !! to the receiver nodes.
    !!
    !! These matrices get stored in:
    !!      mesh(idom)%chimera%recv(ChiID)%donor_interpolator
    !!      mesh(idom)%chimera%recv(ChiID)%donor_interpolator_grad1
    !!      mesh(idom)%chimera%recv(ChiID)%donor_interpolator_grad2
    !!      mesh(idom)%chimera%recv(ChiID)%donor_interpolator_grad3
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-----------------------------------------------------------------------------------------------------------------
    subroutine compute_chimera_interpolators(mesh)
        type(mesh_t),   intent(inout)   :: mesh

        integer(ik)     :: idom, ChiID, idonor, ierr, ipt, iterm,   &
                           donor_idomain_g, donor_idomain_l, donor_ielement_g, donor_ielement_l, &
                           npts, donor_nterms_s, spacedim
        real(rk)        :: node(3)

        real(rk)        :: jinv, ddxi, ddeta, ddzeta
        real(rk), allocatable, dimension(:,:)   ::  &
            interpolator, interpolator_grad1, interpolator_grad2, interpolator_grad3, metric

        

        !
        ! Loop over all domains
        !
        do idom = 1,mesh%ndomains()

            !
            ! Loop over each chimera face
            !
            do ChiID = 1,mesh%domain(idom)%chimera%nreceivers()

                
                !
                ! For each donor, compute an interpolation matrix
                !
                do idonor = 1,mesh%domain(idom)%chimera%recv(ChiID)%ndonors()

                    donor_idomain_g  = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%idomain_g
                    donor_idomain_l  = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%idomain_l
                    donor_ielement_g = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%ielement_g
                    donor_ielement_l = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%ielement_l
                    donor_nterms_s   = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%nterms_s

                    !
                    ! Get number of GQ points this donor is responsible for
                    !
                    npts   = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%nnodes()

                    !
                    ! Allocate interpolator matrix
                    !
                    if (allocated(interpolator)) deallocate(interpolator,       &
                                                            interpolator_grad1, &
                                                            interpolator_grad2, &
                                                            interpolator_grad3)
                    allocate(interpolator(      npts,donor_nterms_s), &
                             interpolator_grad1(npts,donor_nterms_s), &
                             interpolator_grad2(npts,donor_nterms_s), &
                             interpolator_grad3(npts,donor_nterms_s), stat=ierr)
                    if (ierr /= 0) call AllocationError

                    !
                    ! Compute values of modal polynomials at the donor nodes
                    !
                    do iterm = 1,donor_nterms_s
                        do ipt = 1,npts

                            node = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%coords(ipt,:)

                            !
                            ! Compute value interpolator
                            !
                            spacedim = 3
                            interpolator(ipt,iterm) = polynomial_val(spacedim,donor_nterms_s,iterm,node)

                            
                            !
                            ! Compute gradient interpolators, grad1, grad2, grad3
                            !
                            ddxi   = dpolynomial_val(spacedim,donor_nterms_s,iterm,node,XI_DIR  )
                            ddeta  = dpolynomial_val(spacedim,donor_nterms_s,iterm,node,ETA_DIR )
                            ddzeta = dpolynomial_val(spacedim,donor_nterms_s,iterm,node,ZETA_DIR)

                            ! Get metrics for element mapping
                            metric = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%metric(:,:,ipt)
                            jinv   = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%jinv(ipt)

                            ! Compute cartesian derivative interpolator for gq node
                            interpolator_grad1(ipt,iterm) = metric(1,1) * ddxi   + &
                                                            metric(2,1) * ddeta  + &
                                                            metric(3,1) * ddzeta 
                            interpolator_grad2(ipt,iterm) = metric(1,2) * ddxi   + &
                                                            metric(2,2) * ddeta  + &
                                                            metric(3,2) * ddzeta 
                            interpolator_grad3(ipt,iterm) = metric(1,3) * ddxi   + &
                                                            metric(2,3) * ddeta  + &
                                                            metric(3,3) * ddzeta 

                        end do ! ipt
                    end do ! iterm

                    !
                    ! Store interpolators
                    !
                    mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%value = interpolator
                    mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%grad1 = interpolator_grad1
                    mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%grad2 = interpolator_grad2
                    mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%grad3 = interpolator_grad3


                end do  ! idonor



            end do  ! ChiID
        end do  ! idom



        
        !
        ! Communicate mesh
        !
        call mesh%comm_send()
        call mesh%comm_recv()   ! also calls mesh%assemble_chimera_data to construct data on complete exterior node set.
        call mesh%comm_wait()



    end subroutine compute_chimera_interpolators
    !******************************************************************************************************************











    !>  Find the domain and element indices for an element that contains a given quadrature node and can donate
    !!  interpolated solution values to the receiver face.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]      mesh                Array of mesh_t instances
    !!  @param[in]      gq_node             GQ point that needs to find a donor
    !!  @param[in]      receiver_face       Location of face containing the gq_node
    !!  @param[inout]   donor_element       Location of the donor element that was found
    !!  @param[inout]   donor_coordinate    Point defining the location of the GQ point in the donor coordinate system
    !!  @param[inout]   donor_volume        Volume of the donor element that can be used to select between donors if 
    !!                                      multiple are available.
    !!
    !-----------------------------------------------------------------------------------------------------------------------
    !subroutine find_gq_donor(mesh,gq_node,receiver_face,donor_element,donor_coordinate,donor_volume)
    subroutine find_gq_donor(mesh,gq_node,receiver_face,donor_element,donor_coordinate,donor_found,donor_volume)
        type(mesh_t),               intent(in)              :: mesh
        type(point_t),              intent(in)              :: gq_node
        type(face_info_t),          intent(in)              :: receiver_face
        type(element_info_t),       intent(inout)           :: donor_element
        !type(point_t),              intent(inout)           :: donor_coordinate
        real(rk),                   intent(inout)           :: donor_coordinate(3)
        logical,                    intent(inout)           :: donor_found
        real(rk),                   intent(inout), optional :: donor_volume


        integer(ik), allocatable    :: domains_g(:)
        integer(ik)                 :: idom, ielem, inewton, idomain_g, idomain_l,      &
                                       ielement_g, ielement_l, icandidate, ncandidates, &
                                       idonor, ndonors, donor_index

        real(rk), allocatable   :: xcenter(:), ycenter(:), zcenter(:), dist(:), donor_vols(:)
        real(rk)                :: xgq, ygq, zgq, dx, dy, dz, xi, eta, zeta, xn, yn, zn,    &
                                   xmin, xmax, ymin, ymax, zmin, zmax,                      &
                                   xcenter_recv, ycenter_recv, zcenter_recv

        real(rk)                :: gq_comp(3)
        type(ivector_t)         :: candidate_domains_g, candidate_domains_l, candidate_elements_g, candidate_elements_l
        type(ivector_t)         :: donors
        type(rvector_t)         :: donors_xi, donors_eta, donors_zeta

        logical                 :: contained = .false.
        logical                 :: receiver  = .false.
        logical                 :: node_found = .false.



        xgq = gq_node%c1_
        ygq = gq_node%c2_
        zgq = gq_node%c3_


        !
        ! Loop through LOCAL domains and search for potential donor candidates
        !
        ncandidates = 0
        do idom = 1,mesh%ndomains()
            idomain_g = mesh%domain(idom)%idomain_g
            idomain_l = mesh%domain(idom)%idomain_l



            !
            ! Loop through elements in the current domain
            !
            do ielem = 1,mesh%domain(idom)%nelem
                ielement_g = mesh%domain(idom)%elems(ielem)%ielement_g
                ielement_l = mesh%domain(idom)%elems(ielem)%ielement_l

                !
                ! Get bounding coordinates for the current element
                !
                xmin = minval(mesh%domain(idom)%elems(ielem)%elem_pts(:,1))
                xmax = maxval(mesh%domain(idom)%elems(ielem)%elem_pts(:,1))
                ymin = minval(mesh%domain(idom)%elems(ielem)%elem_pts(:,2))
                ymax = maxval(mesh%domain(idom)%elems(ielem)%elem_pts(:,2))
                zmin = minval(mesh%domain(idom)%elems(ielem)%elem_pts(:,3))
                zmax = maxval(mesh%domain(idom)%elems(ielem)%elem_pts(:,3))

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
            gq_comp = mesh%domain(idomain_l)%elems(ielement_l)%computational_point(xgq,ygq,zgq)    ! Newton's method routine
            node_found = (any(ieee_is_nan(gq_comp)) .eqv. .false.)

            ! Add donor if gq_comp point is valid
            !
            ! NOTE: the exit call here enforces that only one candidate is considered and others are thrown away. Maybe not the best choice.
            !
            !donor_found = (gq_comp%status == VALID_POINT)
            !donor_found = (donor_status == VALID_POINT)

            if ( node_found ) then
                ndonors = ndonors + 1
                call donors%push_back(icandidate)
                call donors_xi%push_back(  gq_comp(1))
                call donors_eta%push_back( gq_comp(2))
                call donors_zeta%push_back(gq_comp(3))
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

            !donor_coordinate%status = INVALID_POINT
            donor_found = .false.



        elseif (ndonors == 1) then
            idonor = donors%at(1)   ! donor index from candidates
            donor_element%idomain_g  = candidate_domains_g%at(idonor)
            donor_element%idomain_l  = candidate_domains_l%at(idonor)
            donor_element%ielement_g = candidate_elements_g%at(idonor)
            donor_element%ielement_l = candidate_elements_l%at(idonor)
            donor_element%iproc      = IRANK
            donor_element%eqn_ID     = mesh%domain(donor_element%idomain_l)%eqn_ID
            donor_element%neqns      = mesh%domain(donor_element%idomain_l)%elems(donor_element%ielement_l)%neqns
            donor_element%nterms_s   = mesh%domain(donor_element%idomain_l)%elems(donor_element%ielement_l)%nterms_s
            donor_element%nterms_c   = mesh%domain(donor_element%idomain_l)%elems(donor_element%ielement_l)%nterms_c

            xi   = donors_xi%at(1)
            eta  = donors_eta%at(1)
            zeta = donors_zeta%at(1)
            !call donor_coordinate%set(xi,eta,zeta)
            !donor_coordinate%status = VALID_POINT
            donor_coordinate = [xi,eta,zeta]
            donor_found = .true.
            if (present(donor_volume)) donor_volume = mesh%domain(donor_element%idomain_l)%elems(donor_element%ielement_l)%vol



        elseif (ndonors > 1) then
            !
            ! Handle multiple potential donors: Choose donor with minimum volume - should be best resolved
            !
            if (allocated(donor_vols) ) deallocate(donor_vols)
            allocate(donor_vols(donors%size()))
            
            do idonor = 1,donors%size()
                donor_vols(idonor) = mesh%domain(candidate_domains_l%at(donors%at(idonor)))%elems(candidate_elements_l%at(donors%at(idonor)))%vol
            end do 
    

            !
            ! Get index of domain with minimum volume
            !
            donor_index = minloc(donor_vols,1)
            idonor = donors%at(donor_index)

            donor_element%idomain_g  = candidate_domains_g%at(idonor)
            donor_element%idomain_l  = candidate_domains_l%at(idonor)
            donor_element%ielement_g = candidate_elements_g%at(idonor)
            donor_element%ielement_l = candidate_elements_l%at(idonor)
            donor_element%iproc      = IRANK
            donor_element%eqn_ID     = mesh%domain(donor_element%idomain_l)%eqn_ID
            donor_element%neqns      = mesh%domain(donor_element%idomain_l)%elems(donor_element%ielement_l)%neqns
            donor_element%nterms_s   = mesh%domain(donor_element%idomain_l)%elems(donor_element%ielement_l)%nterms_s
            donor_element%nterms_c   = mesh%domain(donor_element%idomain_l)%elems(donor_element%ielement_l)%nterms_c

            !
            ! Set donor coordinate and volume if present
            !
            xi   = donors_xi%at(donor_index)
            eta  = donors_eta%at(donor_index)
            zeta = donors_zeta%at(donor_index)
            !call donor_coordinate%set(xi,eta,zeta)
            !donor_coordinate%status = VALID_POINT
            donor_coordinate = [xi,eta,zeta]
            donor_found = .true.
            if (present(donor_volume)) donor_volume = mesh%domain(donor_element%idomain_l)%elems(donor_element%ielement_l)%vol


        else
            call chidg_signal(FATAL,"find_gq_donor: invalid number of donors")
        end if




    end subroutine find_gq_donor
    !****************************************************************************************************************









end module mod_chimera
