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
                                      ONE, ZERO, TWO_DIM, THREE_DIM, RKTOL

    use type_mesh,              only: mesh_t
    use type_point,             only: point_t
    use type_element_indices,   only: element_indices_t
    use type_face_info,         only: face_info_t
    use type_ivector,           only: ivector_t
    use type_rvector,           only: rvector_t
    use type_pvector,           only: pvector_t

    use mod_polynomial,         only: polynomialVal
    use mod_periodic,           only: compute_periodic_offset
    use mod_inv,                only: inv
    use mod_chidg_mpi,          only: IRANK, NRANK
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

        integer(ik) :: idom, igq, ichimera_face, idonor, ierr
        integer(ik) :: ndonors, neqns, nterms_s
        integer(ik) :: idonor_domain_g, idonor_element_g
        integer(ik) :: idonor_domain_l, idonor_element_l
        integer(ik) :: idomain_g_list, idomain_l_list, ielement_g_list, ielement_l_list

        real(rk)    :: offset_x, offset_y, offset_z

        type(face_info_t)           :: receiver
        type(element_indices_t)     :: donor
        type(point_t)               :: donor_coord
        type(point_t)               :: gq_node
        type(point_t)               :: dummy_coord
        logical                     :: new_donor     = .false.
        logical                     :: already_added = .false.
        logical                     :: donor_match   = .false.

        type(ivector_t)             :: ddomain_g, ddomain_l, delement_g, delement_l
        type(pvector_t)             :: dcoordinate

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
                receiver%iface    = mesh(idom)%chimera%recv%data(ichimera_face)%receiver_face

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


                    !
                    ! Call routine to find gq donor for current node
                    !
                    call compute_gq_donor(mesh,gq_node, receiver, donor, donor_coord)


                    !
                    ! Add donor location and coordinate
                    !
                    call ddomain_g%push_back(donor%idomain_g)
                    call ddomain_l%push_back(donor%idomain_l)
                    call delement_g%push_back(donor%ielement_g)
                    call delement_l%push_back(donor%ielement_l)
                    call dcoordinate%push_back(donor_coord)


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


                    !
                    ! Check if domain/element pair has already been added to the chimera donor data
                    !
                    already_added = .false.
                    do idonor = 1,mesh(idom)%chimera%recv%data(ichimera_face)%donor_domain_g%size()

                        idonor_domain_g  = mesh(idom)%chimera%recv%data(ichimera_face)%donor_domain_g%at(idonor)
                        idonor_domain_l  = mesh(idom)%chimera%recv%data(ichimera_face)%donor_domain_l%at(idonor)
                        idonor_element_g = mesh(idom)%chimera%recv%data(ichimera_face)%donor_element_g%at(idonor)
                        idonor_element_l = mesh(idom)%chimera%recv%data(ichimera_face)%donor_element_l%at(idonor)

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
                        neqns    = mesh(idomain_l_list)%elems(ielement_l_list)%neqns
                        nterms_s = mesh(idomain_l_list)%elems(ielement_l_list)%nterms_s

                        call mesh(idom)%chimera%recv%data(ichimera_face)%donor_neqns%push_back(neqns)
                        call mesh(idom)%chimera%recv%data(ichimera_face)%donor_nterms_s%push_back(nterms_s)
                        call mesh(idom)%chimera%recv%data(ichimera_face)%donor_domain_g%push_back(idomain_g_list)
                        call mesh(idom)%chimera%recv%data(ichimera_face)%donor_domain_l%push_back(idomain_l_list)
                        call mesh(idom)%chimera%recv%data(ichimera_face)%donor_element_g%push_back(ielement_g_list)
                        call mesh(idom)%chimera%recv%data(ichimera_face)%donor_element_l%push_back(ielement_l_list)
                        ndonors = ndonors + 1
                    end if

                end do ! igq



                !
                ! Allocate chimera donor coordinate and quadrature index arrays. One list for each donor
                !
                mesh(idom)%chimera%recv%data(ichimera_face)%ndonors = ndonors
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
                ! Clear temporary face arrays
                !
                call ddomain_g%clear()
                call ddomain_l%clear()
                call delement_g%clear()
                call delement_l%clear()
                call dcoordinate%clear()


            end do ! iface

        end do ! idom



    end subroutine detect_chimera_donors
    !***********************************************************************************************************************









    !> compute the donor domain and element for a given quadrature node.
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
    subroutine compute_gq_donor(mesh,gq_node,receiver_face,donor_element,donor_coordinate)
        type(mesh_t),               intent(in)      :: mesh(:)
        type(point_t),              intent(in)      :: gq_node
        type(face_info_t),          intent(in)      :: receiver_face
        type(element_indices_t),    intent(inout)   :: donor_element
        type(point_t),              intent(inout)   :: donor_coordinate


        integer(ik)             :: idom, ielem, inewton, spacedim
        integer(ik)             :: idomain_g, idomain_l, ielement_g, ielement_l
        integer(ik)             :: icandidate, ncandidates, idonor, ndonors
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
        logical                 :: contained = .false.
        logical                 :: receiver  = .false.
        logical                 :: valid_point



        tol = 10._rk*RKTOL


        xgq = gq_node%c1_
        ygq = gq_node%c2_
        zgq = gq_node%c3_


        !
        ! Loop through domains and search for potential donor candidates
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
                receiver = ( (idomain_g == receiver_face%idomain_g) .and. (ielement_g == receiver_face%ielement_g) .and. &
                             (idomain_l == receiver_face%idomain_l) .and. (ielement_l == receiver_face%ielement_l) )

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
        ! Test gq_node on candidate element volume using Newton's method to map to donor local coordinates
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
            gq_comp = mesh(idomain_l)%elems(ielement_l)%computational_point(xgq,ygq,zgq)

            !
            ! Add donor if gq_comp point is valid
            !
            valid_point = (gq_comp%status == 0)
            if ( valid_point ) then
                ndonors = ndonors + 1
                call donors%push_back(icandidate)
                call donors_xi%push_back(  gq_comp%c1_)
                call donors_eta%push_back( gq_comp%c2_)
                call donors_zeta%push_back(gq_comp%c3_)
                exit
            end if

        end do ! icandidate




        !
        ! Sanity check on donors and set donor_element location
        !
        if (ndonors == 0) then
            call chidg_signal_three(FATAL,"compute_gq_donor: No donor found for gq_node", gq_node%c1_, gq_node%c2_, gq_node%c3_)

        elseif (ndonors > 1) then
            !TODO: Account for case of multiple overlapping donors. When a gq node could be filled by two or more elements.
            !      Maybe, just choose one. Maybe, average contribution from all potential donors.
            !
            call chidg_signal(FATAL,"compute_gq_donor: Multiple donors found for the same gq_node")

        elseif (ndonors == 1) then
            idonor = donors%at(1)   ! donor index from candidates
            donor_element%idomain_g  = candidate_domains_g%at(idonor)
            donor_element%idomain_l  = candidate_domains_l%at(idonor)
            donor_element%ielement_g = candidate_elements_g%at(idonor)
            donor_element%ielement_l = candidate_elements_l%at(idonor)

            xi   = donors_xi%at(1)
            eta  = donors_eta%at(1)
            zeta = donors_zeta%at(1)
            call donor_coordinate%set(xi,eta,zeta)

        else
            call chidg_signal(FATAL,"compute_gq_donor: invalid number of donors")
        end if




    end subroutine compute_gq_donor
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
        integer(ik) :: npts, nterms_s, nterms, spacedim

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
                do idonor = 1,mesh(idom)%chimera%recv%data(iChiID)%ndonors

                    donor_idomain_g  = mesh(idom)%chimera%recv%data(iChiID)%donor_domain_g%at(idonor)
                    donor_idomain_l  = mesh(idom)%chimera%recv%data(iChiID)%donor_domain_l%at(idonor)
                    donor_ielement_g = mesh(idom)%chimera%recv%data(iChiID)%donor_element_g%at(idonor)
                    donor_ielement_l = mesh(idom)%chimera%recv%data(iChiID)%donor_element_l%at(idonor)

                    !
                    ! Get number of GQ points this donor is responsible for
                    !
                    npts   = mesh(idom)%chimera%recv%data(iChiID)%donor_coords(idonor)%size()
                    nterms = mesh(donor_idomain_l)%elems(donor_ielement_l)%nterms_s

                    !
                    ! Allocate interpolator matrix
                    !
                    if (allocated(interpolator)) deallocate(interpolator)
                    allocate(interpolator(npts,nterms), stat=ierr)
                    if (ierr /= 0) call AllocationError

                    !
                    ! Compute values of modal polynomials at the donor nodes
                    !
                    do iterm = 1,nterms
                        do ipt = 1,npts

                            node = mesh(idom)%chimera%recv%data(iChiID)%donor_coords(idonor)%at(ipt)
                            interpolator(ipt,iterm) = polynomialVal(spacedim,nterms,iterm,node)

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
