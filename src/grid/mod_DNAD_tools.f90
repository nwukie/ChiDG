module mod_DNAD_tools
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, DIAG, CHIMERA, INTERIOR, BOUNDARY
    use mod_chidg_mpi,          only: IRANK
    
    use type_mesh_new,          only: mesh_new_t
    use type_seed,              only: seed_t
    implicit none



contains










    !> Computes the domain index of the neighbor domain
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh    Array of mesh_t instances
    !!  @param[in]  idom    Domain index of the current element
    !!  @param[in]  ielem   Element index for mesh(idom)%elems(ielem)
    !!  @param[in]  iface   Face index for mesh(idom)%faces(ielem,iface)
    !!  @param[in]  idonor  Index of potential donor elements. For example, Chimera faces could have more than one donor element
    !!
    !-----------------------------------------------------------------------------------------------
    function compute_neighbor_domain_l(mesh,idom,ielem,iface,idonor) result(idom_n)
        type(mesh_new_t),   intent(in)  :: mesh
        integer(ik),        intent(in)  :: idom
        integer(ik),        intent(in)  :: ielem
        integer(ik),        intent(in)  :: iface
        integer(ik),        intent(in)  :: idonor

        integer(ik) :: idom_n, ChiID
        logical     :: chimera_face = .false.

        chimera_face = (mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA)


        if ( chimera_face ) then

            ChiID  = mesh%domain(idom)%faces(ielem,iface)%ChiID
            idom_n = mesh%domain(idom)%chimera%recv%data(ChiID)%donor_domain_l%at(idonor)

        else
            idom_n = mesh%domain(idom)%faces(ielem,iface)%ineighbor_domain_l
        end if

    end function compute_neighbor_domain_l
    !************************************************************************************************













    !> Computes the domain index of the neighbor domain
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh    Array of mesh_t instances
    !!  @param[in]  idom    Domain index of the current element
    !!  @param[in]  ielem   Element index for mesh(idom)%elems(ielem)
    !!  @param[in]  iface   Face index for mesh(idom)%faces(ielem,iface)
    !!  @param[in]  idonor  Index of potential donor elements. For example, Chimera faces could have more than one donor element
    !!
    !-----------------------------------------------------------------------------------------------
    function compute_neighbor_domain_g(mesh,idom,ielem,iface,idonor) result(idom_n)
        type(mesh_new_t),   intent(in)  :: mesh
        integer(ik),        intent(in)  :: idom
        integer(ik),        intent(in)  :: ielem
        integer(ik),        intent(in)  :: iface
        integer(ik),        intent(in)  :: idonor

        integer(ik) :: idom_n, ChiID
        logical     :: chimera_face = .false.

        chimera_face = (mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA)


        if ( chimera_face ) then

            ChiID  = mesh%domain(idom)%faces(ielem,iface)%ChiID
            idom_n = mesh%domain(idom)%chimera%recv%data(ChiID)%donor_domain_g%at(idonor)

        else
            idom_n = mesh%domain(idom)%faces(ielem,iface)%ineighbor_domain_g
        end if

    end function compute_neighbor_domain_g
    !************************************************************************************************













    !> Computes the element index of the neighbor element
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh    Array of mesh_t instances
    !!  @param[in]  idom    Domain index of the current element
    !!  @param[in]  ielem   Element index for mesh(idom)%elems(ielem)
    !!  @param[in]  iface   Face index for mesh(idom)%faces(ielem,iface)
    !!  @param[in]  idonor  Index of potential donor elements. For example, Chimera faces could have more than one donor element
    !!
    !------------------------------------------------------------------------------------------------
    function compute_neighbor_element_l(mesh,idom,ielem,iface,idonor) result(ielem_n)
        type(mesh_new_t),   intent(in)  :: mesh
        integer(ik),        intent(in)  :: idom
        integer(ik),        intent(in)  :: ielem
        integer(ik),        intent(in)  :: iface
        integer(ik),        intent(in)  :: idonor

        integer(ik) :: ielem_n, ChiID
        logical     :: chimera_face = .false.

        chimera_face = (mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA)
        
        if ( chimera_face ) then

            ChiID   = mesh%domain(idom)%faces(ielem,iface)%ChiID
            ielem_n = mesh%domain(idom)%chimera%recv%data(ChiID)%donor_element_l%at(idonor)

        else
            ielem_n = mesh%domain(idom)%faces(ielem,iface)%get_neighbor_element_l()
        end if


    end function compute_neighbor_element_l
    !************************************************************************************************








    !> Computes the element index of the neighbor element
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh    Array of mesh_t instances
    !!  @param[in]  idom    Domain index of the current element
    !!  @param[in]  ielem   Element index for mesh(idom)%elems(ielem)
    !!  @param[in]  iface   Face index for mesh(idom)%faces(ielem,iface)
    !!  @param[in]  idonor  Index of potential donor elements. For example, Chimera faces could have more than one donor element
    !!
    !------------------------------------------------------------------------------------------------
    function compute_neighbor_element_g(mesh,idom,ielem,iface,idonor) result(ielem_n)
        type(mesh_new_t),   intent(in)  :: mesh
        integer(ik),        intent(in)  :: idom
        integer(ik),        intent(in)  :: ielem
        integer(ik),        intent(in)  :: iface
        integer(ik),        intent(in)  :: idonor

        integer(ik) :: ielem_n, ChiID
        logical     :: chimera_face = .false.

        chimera_face = (mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA)
        
        if ( chimera_face ) then

            ChiID   = mesh%domain(idom)%faces(ielem,iface)%ChiID
            ielem_n = mesh%domain(idom)%chimera%recv%data(ChiID)%donor_element_l%at(idonor)

        else
            ielem_n = mesh%domain(idom)%faces(ielem,iface)%ineighbor_element_g
        end if


    end function compute_neighbor_element_g
    !************************************************************************************************











    !> Computes the face index of matching face in neighboring element
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh    Array of mesh_t instances
    !!  @param[in]  idom    Domain index of the current element
    !!  @param[in]  ielem   Element index for mesh(idom)%elems(ielem)
    !!  @param[in]  iface   Face index for mesh(idom)%faces(ielem,iface)
    !!  @param[in]  idonor  Index of potential donor elements. For example, Chimera faces could have more than one donor element
    !!
    !---------------------------------------------------------------------------------------------------------
    function compute_neighbor_face(mesh,idom,ielem,iface,idonor) result(iface_n)
        type(mesh_new_t),   intent(in)  :: mesh
        integer(ik),        intent(in)  :: idom
        integer(ik),        intent(in)  :: ielem
        integer(ik),        intent(in)  :: iface
        integer(ik),        intent(in)  :: idonor

        integer(ik) :: iface_n
        logical     :: chimera_face

        chimera_face = ( mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA )

        if ( chimera_face ) then
            iface_n = CHIMERA
        else

            if (iface == XI_MIN) then
                iface_n = XI_MAX
            else if (iface == XI_MAX) then
                iface_n = XI_MIN
            else if (iface == ETA_MIN) then
                iface_n = ETA_MAX
            else if (iface == ETA_MAX) then
                iface_n = ETA_MIN
            else if (iface == ZETA_MIN) then
                iface_n = ZETA_MAX
            else if (iface == ZETA_MAX) then
                iface_n = ZETA_MIN
            end if

        end if

    end function compute_neighbor_face
    !**********************************************************************************************************












    !> Compute the domain and element that is being linearized. These are found by checking
    !! the current element/face neighbors, or by setting the current element itself. Which
    !! element gets linearized depends on idiff. idiff specifies the direction of the linearization.
    !! For example, if idiff == XI_MIN, the seed element is that which is the neighbor to the XI_MIN
    !! face of the current element. The other behavior is if the linearization of the current
    !! element is desired. That is handled with idiff == DIAG.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh    Array of mesh_t instances
    !!  @param[in]  idom    Domain index of the current element
    !!  @param[in]  ielem   Element index for mesh(idom)%elems(ielem)
    !!  @param[in]  iface   Face index for mesh(idom)%faces(ielem,iface)
    !!  @param[in]  idepend Index of external dependent element. For example, Chimera faces could have more than one dependent element
    !!  @param[in]  idiff   Linearization index
    !!
    !-------------------------------------------------------------------------------------------------------------------------
    function face_compute_seed(mesh,idomain_l,ielement_l,iface,idepend,idiff) result(seed)
        type(mesh_new_t),   intent(in)  :: mesh
        integer(ik),        intent(in)  :: idomain_l
        integer(ik),        intent(in)  :: ielement_l
        integer(ik),        intent(in)  :: iface
        integer(ik),        intent(in)  :: idepend
        integer(ik),        intent(in)  :: idiff


        type(seed_t)    :: seed
        integer(ik)     :: ChiID
        logical         :: linearize_me, linearize_neighbor, linearize
        logical         :: chimera_face, interior_face, boundary_face


        linearize           = ( idiff /= 0 )
        linearize_me        = ( idiff == DIAG )
        linearize_neighbor  = ( idiff == XI_MIN   .or. idiff == XI_MAX   .or. &
                                idiff == ETA_MIN  .or. idiff == ETA_MAX  .or. &
                                idiff == ZETA_MIN .or. idiff == ZETA_MAX )


        if (linearize) then

            !
            ! Check for linearization of current element (ielem)
            !
            if ( linearize_me ) then

                call seed%init(idomain_g    = mesh%domain(idomain_l)%elems(ielement_l)%idomain_g,  &
                               idomain_l    = mesh%domain(idomain_l)%elems(ielement_l)%idomain_l,  &
                               ielement_g   = mesh%domain(idomain_l)%elems(ielement_l)%ielement_g, &
                               ielement_l   = mesh%domain(idomain_l)%elems(ielement_l)%ielement_l, &
                               neqns        = mesh%domain(idomain_l)%elems(ielement_l)%neqns,      &
                               nterms_s     = mesh%domain(idomain_l)%elems(ielement_l)%nterms_s,   &
                               iproc        = IRANK,                                        &
                               recv_comm    = 0,                                            &
                               recv_domain  = 0,                                            &
                               recv_element = 0)


            !
            ! Check for linearization of neighbor element (neighbor of face(ielem,iface) )
            !
            elseif ( linearize_neighbor ) then



                !
                ! Check if linearization direction (iface) is a Interior or Chimera face
                !
                chimera_face  = ( mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == CHIMERA  )
                interior_face = ( mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == INTERIOR )
                boundary_face = ( mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == BOUNDARY )


                !
                ! Linearize wrt Standard Interior Neighbor
                !
                if ( interior_face ) then

                    call seed%init(idomain_g    = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_g,   &
                                   idomain_l    = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_l,   &
                                   ielement_g   = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_element_g,  &
                                   ielement_l   = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_element_l,  &
                                   neqns        = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_neqns,      &
                                   nterms_s     = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_nterms_s,   &
                                   iproc        = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_proc,       &
                                   recv_comm    = mesh%domain(idomain_l)%faces(ielement_l,iface)%recv_comm,            &
                                   recv_domain  = mesh%domain(idomain_l)%faces(ielement_l,iface)%recv_domain,          &
                                   recv_element = mesh%domain(idomain_l)%faces(ielement_l,iface)%recv_element)


                !
                ! Linearize wrt Chimera Interior Neighbor
                !
                elseif ( chimera_face ) then
                    ChiID = mesh%domain(idomain_l)%faces(ielement_l,iface)%ChiID

                    call seed%init(idomain_g    = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_domain_g%at(idepend),       &
                                   idomain_l    = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_domain_l%at(idepend),       &
                                   ielement_g   = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_element_g%at(idepend),      &
                                   ielement_l   = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_element_l%at(idepend),      &
                                   neqns        = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_neqns%at(idepend),          &
                                   nterms_s     = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_nterms_s%at(idepend),       &
                                   iproc        = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_proc%at(idepend),           &
                                   recv_comm    = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_recv_comm%at(idepend),      &
                                   recv_domain  = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_recv_domain%at(idepend),    &
                                   recv_element = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_recv_element%at(idepend) )



                !
                ! Boudnary face, linearize wrt Current element
                !
                elseif ( boundary_face ) then

!                    call seed%init(idomain_g    = mesh%domain(idomain_l)%elems(ielement_l)%idomain_g,  &
!                                   idomain_l    = mesh%domain(idomain_l)%elems(ielement_l)%idomain_l,  &
!                                   ielement_g   = mesh%domain(idomain_l)%elems(ielement_l)%ielement_g, &
!                                   ielement_l   = mesh%domain(idomain_l)%elems(ielement_l)%ielement_l, &
!                                   neqns        = mesh%domain(idomain_l)%elems(ielement_l)%neqns,      &
!                                   nterms_s     = mesh%domain(idomain_l)%elems(ielement_l)%nterms_s,   &
!                                   iproc        = IRANK,                                        &
!                                   recv_comm    = 0,                                            &
!                                   recv_domain  = 0,                                            &
!                                   recv_element = 0)


                ! Invalid Case
                else
                    call chidg_signal(FATAL,"face_compute_seed: invalid face type - face(ielem,iface)%ftype")
                end if


            ! Invalid value for iface
            else
                call chidg_signal(FATAL,"face_compute_seed: invalid value for iface")
            end if ! linearize_me or linearize_neighbor


        else

            ! If idiff == 0 then no linearization. 
            call seed%init(idomain_g    = 0, &
                           idomain_l    = 0, &
                           ielement_g   = 0, &
                           ielement_l   = 0, &
                           neqns        = 0, &
                           nterms_s     = 0, &
                           iproc        = 0, &
                           recv_comm    = 0, &
                           recv_domain  = 0, &
                           recv_element = 0)



        end if ! linearize



    end function face_compute_seed
    !************************************************************************************************************












    !> Compute the domain and element that is being linearized. These are found by checking
    !! the current element/face neighbors, or by setting the current element itself. Which
    !! element gets linearized depends on idiff. idiff specifies the direction of the linearization.
    !! For example, if idiff == XI_MIN, the seed element is that which is the neighbor to the XI_MIN
    !! face of the current element. The other behavior is if the linearization of the current
    !! element is desired. That is handled with idiff == DIAG.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh    Array of mesh_t instances
    !!  @param[in]  idom    Domain index of the current element
    !!  @param[in]  ielem   Element index for mesh(idom)%elems(ielem)
    !!  @param[in]  iface   Face index for mesh(idom)%faces(ielem,iface)
    !!  @param[in]  idepend Index of potential dependent exterior elements. For example, Chimera faces could 
    !!                      have more than one dependent element
    !!  @param[in]  idiff   Linearization index
    !!
    !-----------------------------------------------------------------------------------------------------------
    function element_compute_seed(mesh,idomain_l,ielement_l,idepend,idiff) result(seed)
        type(mesh_new_t),   intent(in)  :: mesh
        integer(ik),        intent(in)  :: idomain_l
        integer(ik),        intent(in)  :: ielement_l
        integer(ik),        intent(in)  :: idepend
        integer(ik),        intent(in)  :: idiff


        type(seed_t)    :: seed
        integer(ik)     :: ChiID, iface
        logical         :: linearize_me, linearize_neighbor, linearize
        logical         :: chimera_face, interior_face, boundary_face


        !linearize           = (idepend /= 0)
        linearize           = ( idiff /= 0 )
        linearize_me        = ( idiff == DIAG )
        linearize_neighbor  = ( idiff == XI_MIN   .or. idiff == XI_MAX   .or. &
                                idiff == ETA_MIN  .or. idiff == ETA_MAX  .or. &
                                idiff == ZETA_MIN .or. idiff == ZETA_MAX )


        if (linearize) then

            !
            ! Check for linearization of current element (ielem)
            !
            if ( linearize_me ) then

                call seed%init(idomain_g    = mesh%domain(idomain_l)%elems(ielement_l)%idomain_g,  &
                               idomain_l    = mesh%domain(idomain_l)%elems(ielement_l)%idomain_l,  &
                               ielement_g   = mesh%domain(idomain_l)%elems(ielement_l)%ielement_g, &
                               ielement_l   = mesh%domain(idomain_l)%elems(ielement_l)%ielement_l, &
                               neqns        = mesh%domain(idomain_l)%elems(ielement_l)%neqns,      &
                               nterms_s     = mesh%domain(idomain_l)%elems(ielement_l)%nterms_s,   &
                               iproc        = IRANK,                                        &
                               recv_comm    = 0,                                            &
                               recv_domain  = 0,                                            &
                               recv_element = 0)


            !
            ! Check for linearization of neighbor element (neighbor of face(ielem,iface) )
            !
            elseif ( linearize_neighbor ) then

                !
                ! Check face in the direction of idiff
                !
                iface = idiff

                !
                ! Check if linearization direction (iface) is a Interior or Chimera face
                !
                chimera_face  = ( mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == CHIMERA  )
                interior_face = ( mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == INTERIOR )
                boundary_face = ( mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == BOUNDARY )


                !
                ! Linearize wrt Standard Interior Neighbor
                !
                if ( interior_face ) then

                    call seed%init(idomain_g    = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_g,   &
                                   idomain_l    = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_l,   &
                                   ielement_g   = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_element_g,  &
                                   ielement_l   = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_element_l,  &
                                   neqns        = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_neqns,      &
                                   nterms_s     = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_nterms_s,   &
                                   iproc        = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_proc,       &
                                   recv_comm    = mesh%domain(idomain_l)%faces(ielement_l,iface)%recv_comm,            &
                                   recv_domain  = mesh%domain(idomain_l)%faces(ielement_l,iface)%recv_domain,          &
                                   recv_element = mesh%domain(idomain_l)%faces(ielement_l,iface)%recv_element)


                !
                ! Linearize wrt Chimera Interior Neighbor
                !
                elseif ( chimera_face ) then
                    ChiID = mesh%domain(idomain_l)%faces(ielement_l,iface)%ChiID

                    call seed%init(idomain_g    = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_domain_g%at(idepend),       &
                                   idomain_l    = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_domain_l%at(idepend),       &
                                   ielement_g   = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_element_g%at(idepend),      &
                                   ielement_l   = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_element_l%at(idepend),      &
                                   neqns        = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_neqns%at(idepend),          &
                                   nterms_s     = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_nterms_s%at(idepend),       &
                                   iproc        = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_proc%at(idepend),           &
                                   recv_comm    = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_recv_comm%at(idepend),      &
                                   recv_domain  = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_recv_domain%at(idepend),    &
                                   recv_element = mesh%domain(idomain_l)%chimera%recv%data(ChiID)%donor_recv_element%at(idepend) )


                !
                ! Boudnary face, linearize wrt Current element
                !
                elseif ( boundary_face ) then

                    call seed%init(idomain_g    = mesh%domain(idomain_l)%elems(ielement_l)%idomain_g,  &
                                   idomain_l    = mesh%domain(idomain_l)%elems(ielement_l)%idomain_l,  &
                                   ielement_g   = mesh%domain(idomain_l)%elems(ielement_l)%ielement_g, &
                                   ielement_l   = mesh%domain(idomain_l)%elems(ielement_l)%ielement_l, &
                                   neqns        = mesh%domain(idomain_l)%elems(ielement_l)%neqns,      &
                                   nterms_s     = mesh%domain(idomain_l)%elems(ielement_l)%nterms_s,   &
                                   iproc        = IRANK,                                        &
                                   recv_comm    = 0,                                            &
                                   recv_domain  = 0,                                            &
                                   recv_element = 0)


                !
                ! Invalid Case
                !
                else
                    call chidg_signal(FATAL,"element_compute_seed: invalid face type - face(ielem,iface)%ftype")
                end if


            else
                call chidg_signal(FATAL,"element_compute_seed: invalid value for iface")
            end if

        else

            ! If idiff == 0 then no linearization. 
            call seed%init(idomain_g    = 0, &
                           idomain_l    = 0, &
                           ielement_g   = 0, &
                           ielement_l   = 0, &
                           neqns        = 0, &
                           nterms_s     = 0, &
                           iproc        = 0, &
                           recv_comm    = 0, &
                           recv_domain  = 0, &
                           recv_element = 0)


        end if




    end function element_compute_seed
    !************************************************************************************************************























end module mod_DNAD_tools
