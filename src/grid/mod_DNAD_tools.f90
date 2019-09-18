module mod_DNAD_tools
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, DIAG, CHIMERA, INTERIOR, BOUNDARY
    use mod_chidg_mpi,          only: IRANK
    
    use type_mesh,          only: mesh_t
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
        type(mesh_t),   intent(in)  :: mesh
        integer(ik),        intent(in)  :: idom
        integer(ik),        intent(in)  :: ielem
        integer(ik),        intent(in)  :: iface
        integer(ik),        intent(in)  :: idonor

        integer(ik) :: idom_n, ChiID
        logical     :: chimera_face = .false.

        chimera_face = (mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA)

        if ( chimera_face ) then
            ChiID  = mesh%domain(idom)%faces(ielem,iface)%ChiID
            idom_n = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%elem_info%idomain_l
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
        type(mesh_t),   intent(in)  :: mesh
        integer(ik),        intent(in)  :: idom
        integer(ik),        intent(in)  :: ielem
        integer(ik),        intent(in)  :: iface
        integer(ik),        intent(in)  :: idonor

        integer(ik) :: idom_n, ChiID
        logical     :: chimera_face = .false.

        chimera_face = (mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA)

        if ( chimera_face ) then
            ChiID  = mesh%domain(idom)%faces(ielem,iface)%ChiID
            idom_n = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%elem_info%idomain_g
        else
            idom_n = mesh%domain(idom)%faces(ielem,iface)%ineighbor_domain_g
        end if

    end function compute_neighbor_domain_g
    !************************************************************************************************








    !> Computes the local element index of the neighbor element
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
        type(mesh_t),   intent(in)  :: mesh
        integer(ik),        intent(in)  :: idom
        integer(ik),        intent(in)  :: ielem
        integer(ik),        intent(in)  :: iface
        integer(ik),        intent(in)  :: idonor

        integer(ik) :: ielem_n, ChiID
        logical     :: chimera_face = .false.

        chimera_face = (mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA)
        
        if ( chimera_face ) then
            ChiID   = mesh%domain(idom)%faces(ielem,iface)%ChiID
            ielem_n = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%elem_info%ielement_l
        else
            ielem_n = mesh%domain(idom)%faces(ielem,iface)%get_neighbor_element_l()
        end if

    end function compute_neighbor_element_l
    !************************************************************************************************





    !> Computes the global element index of the neighbor element
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
        type(mesh_t),   intent(in)  :: mesh
        integer(ik),        intent(in)  :: idom
        integer(ik),        intent(in)  :: ielem
        integer(ik),        intent(in)  :: iface
        integer(ik),        intent(in)  :: idonor

        integer(ik) :: ielem_n, ChiID
        logical     :: chimera_face = .false.

        chimera_face = (mesh%domain(idom)%faces(ielem,iface)%ftype == CHIMERA)
        
        if ( chimera_face ) then
            ChiID   = mesh%domain(idom)%faces(ielem,iface)%ChiID
            ielem_n = mesh%domain(idom)%chimera%recv(ChiID)%donor(idonor)%elem_info%ielement_g
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
        type(mesh_t),   intent(in)  :: mesh
        integer(ik),    intent(in)  :: idom
        integer(ik),    intent(in)  :: ielem
        integer(ik),    intent(in)  :: iface
        integer(ik),    intent(in)  :: idonor

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
    function face_compute_seed(mesh,idomain_l,ielement_l,iface,idepend,idiff,itime) result(seed)
        type(mesh_t),   intent(in)  :: mesh
        integer(ik),    intent(in)  :: idomain_l
        integer(ik),    intent(in)  :: ielement_l
        integer(ik),    intent(in)  :: iface
        integer(ik),    intent(in)  :: idepend
        integer(ik),    intent(in)  :: idiff
        integer(ik),    intent(in)  :: itime


        type(seed_t)    :: seed
        integer(ik)     :: ChiID, group_ID, patch_ID, face_ID
        logical         :: linearize_interior, linearize_exterior, linearize
        logical         :: chimera_face, interior_face, boundary_face


        linearize          = ( idiff /= 0 )
        linearize_interior = ( idiff == DIAG )
        linearize_exterior = ( idiff == XI_MIN   .or. idiff == XI_MAX   .or. &
                               idiff == ETA_MIN  .or. idiff == ETA_MAX  .or. &
                               idiff == ZETA_MIN .or. idiff == ZETA_MAX )

        if (linearize) then

            ! Check for linearization of current element (ielem)
            if ( linearize_interior ) then

                call seed%init(idomain_g    = mesh%domain(idomain_l)%elems(ielement_l)%idomain_g,   &
                               idomain_l    = mesh%domain(idomain_l)%elems(ielement_l)%idomain_l,   &
                               ielement_g   = mesh%domain(idomain_l)%elems(ielement_l)%ielement_g,  &
                               ielement_l   = mesh%domain(idomain_l)%elems(ielement_l)%ielement_l,  &
                               nfields      = mesh%domain(idomain_l)%elems(ielement_l)%nfields,     &
                               nterms_s     = mesh%domain(idomain_l)%elems(ielement_l)%nterms_s,    &
                               iproc        = IRANK,                                                &
                               itime        = itime,                                                &
                               dof_start    = mesh%domain(idomain_l)%elems(ielement_l)%dof_start,   &
                               recv_comm    = 0,                                                    &
                               recv_domain  = 0,                                                    &
                               recv_element = 0)


            ! Check for linearization of neighbor element (neighbor of face(ielem,iface) )
            elseif ( linearize_exterior ) then

                ! Check if linearization direction (iface) is a Interior or Chimera face
                chimera_face  = ( mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == CHIMERA  )
                interior_face = ( mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == INTERIOR )
                boundary_face = ( mesh%domain(idomain_l)%faces(ielement_l,iface)%ftype == BOUNDARY )

                ! Linearize wrt Standard Interior Neighbor
                if ( interior_face ) then

                    call seed%init(idomain_g    = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_g,    &
                                   idomain_l    = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_l,    &
                                   ielement_g   = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_element_g,   &
                                   ielement_l   = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_element_l,   &
                                   nfields      = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_nfields,     &
                                   nterms_s     = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_nterms_s,    &
                                   iproc        = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_proc,        &
                                   itime        = itime,                                                                &
                                   dof_start    = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_dof_start,   &
                                   recv_comm    = mesh%domain(idomain_l)%faces(ielement_l,iface)%recv_comm,             &
                                   recv_domain  = mesh%domain(idomain_l)%faces(ielement_l,iface)%recv_domain,           &
                                   recv_element = mesh%domain(idomain_l)%faces(ielement_l,iface)%recv_element)

                ! Linearize wrt Chimera Interior Neighbor
                elseif ( chimera_face ) then
                    ChiID = mesh%domain(idomain_l)%faces(ielement_l,iface)%ChiID

                    call seed%init(idomain_g    = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%idomain_g,      &
                                   idomain_l    = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%idomain_l,      &
                                   ielement_g   = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%ielement_g,     &
                                   ielement_l   = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%ielement_l,     &
                                   nfields      = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%nfields,        &
                                   nterms_s     = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%nterms_s,       &
                                   iproc        = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%iproc,          &
                                   itime        = itime,                                                                              &
                                   dof_start    = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%dof_start,      &
                                   recv_comm    = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%recv_comm,      &
                                   recv_domain  = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%recv_domain,    &
                                   recv_element = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%recv_element )



                !
                ! Boudnary face, linearize wrt Current element
                !
                elseif ( boundary_face ) then

                    group_ID = mesh%domain(idomain_l)%faces(ielement_l,iface)%group_ID
                    patch_ID = mesh%domain(idomain_l)%faces(ielement_l,iface)%patch_ID
                    face_ID  = mesh%domain(idomain_l)%faces(ielement_l,iface)%face_ID

                    call seed%init(idomain_g    = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_g( idepend),  &
                                   idomain_l    = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_l( idepend),  &
                                   ielement_g   = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_g(idepend),  &
                                   ielement_l   = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_l(idepend),  &
                                   nfields      = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%nfields(idepend),     &
                                   nterms_s     = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%nterms_s(idepend),    &
                                   iproc        = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%proc(idepend),        &
                                   itime        = itime,                                                                                &
                                   dof_start    = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%dof_start(idepend),   &
                                   recv_comm    = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%recv_comm(idepend),   &
                                   recv_domain  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%recv_domain(idepend), &
                                   recv_element = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%recv_element(idepend) )



                ! Invalid Case
                else
                    call chidg_signal(FATAL,"face_compute_seed: invalid face type - face(ielem,iface)%ftype")
                end if


            ! Invalid value for iface
            else
                call chidg_signal(FATAL,"face_compute_seed: invalid value for iface")
            end if ! linearize_interior or linearize_exterior


        else

            ! no linearization. 
            call seed%init(idomain_g    = 0,            &
                           idomain_l    = 0,            &
                           ielement_g   = 0,            &
                           ielement_l   = 0,            &
                           nfields      = 0,            &
                           nterms_s     = 0,            &
                           iproc        = 0,            &
                           itime        = itime,        &  ! need itime here because the residual storage relies on it
                           dof_start    = 0,            &
                           recv_comm    = 0,            &
                           recv_domain  = 0,            &
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
    function element_compute_seed(mesh,idomain_l,ielement_l,idepend,idiff,itime_couple) result(seed)
        type(mesh_t),   intent(in)  :: mesh
        integer(ik),        intent(in)  :: idomain_l
        integer(ik),        intent(in)  :: ielement_l
        integer(ik),        intent(in)  :: idepend
        integer(ik),        intent(in)  :: idiff
        integer(ik),        intent(in)  :: itime_couple


        type(seed_t)    :: seed
        integer(ik)     :: ChiID, group_ID, patch_ID, face_ID, iface
        logical         :: linearize_interior, linearize_exterior, linearize
        logical         :: chimera_face, interior_face, boundary_face


        linearize           = ( idiff /= 0 )
        linearize_interior  = ( idiff == DIAG )
        linearize_exterior  = ( idiff == XI_MIN   .or. idiff == XI_MAX   .or. &
                                idiff == ETA_MIN  .or. idiff == ETA_MAX  .or. &
                                idiff == ZETA_MIN .or. idiff == ZETA_MAX )


        if (linearize) then

            !
            ! Check for linearization of current element (ielem)
            !
            if ( linearize_interior ) then

                call seed%init(idomain_g    = mesh%domain(idomain_l)%elems(ielement_l)%idomain_g,   &
                               idomain_l    = mesh%domain(idomain_l)%elems(ielement_l)%idomain_l,   &
                               ielement_g   = mesh%domain(idomain_l)%elems(ielement_l)%ielement_g,  &
                               ielement_l   = mesh%domain(idomain_l)%elems(ielement_l)%ielement_l,  &
                               nfields      = mesh%domain(idomain_l)%elems(ielement_l)%nfields,     &
                               nterms_s     = mesh%domain(idomain_l)%elems(ielement_l)%nterms_s,    &
                               iproc        = IRANK,                                                &
                               itime        = itime_couple,                                         &
                               dof_start    = mesh%get_dof_start(idomain_l,ielement_l),             &
                               recv_comm    = 0,                                                    &
                               recv_domain  = 0,                                                    &
                               recv_element = 0)


            !
            ! Check for linearization of exterior element (neighbor of face(ielem,iface) )
            !
            elseif ( linearize_exterior ) then

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

                    call seed%init(idomain_g    = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_g,    &
                                   idomain_l    = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_l,    &
                                   ielement_g   = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_element_g,   &
                                   ielement_l   = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_element_l,   &
                                   nfields      = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_nfields,     &
                                   nterms_s     = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_nterms_s,    &
                                   iproc        = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_proc,        &
                                   itime        = itime_couple,                                                         &
                                   dof_start    = mesh%domain(idomain_l)%faces(ielement_l,iface)%ineighbor_dof_start,   &
                                   recv_comm    = mesh%domain(idomain_l)%faces(ielement_l,iface)%recv_comm,             &
                                   recv_domain  = mesh%domain(idomain_l)%faces(ielement_l,iface)%recv_domain,           &
                                   recv_element = mesh%domain(idomain_l)%faces(ielement_l,iface)%recv_element)


                !
                ! Linearize wrt Chimera Interior Neighbor
                !
                elseif ( chimera_face ) then
                    ChiID = mesh%domain(idomain_l)%faces(ielement_l,iface)%ChiID

                    call seed%init(idomain_g    = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%idomain_g,      &
                                   idomain_l    = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%idomain_l,      &
                                   ielement_g   = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%ielement_g,     &
                                   ielement_l   = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%ielement_l,     &
                                   nfields      = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%nfields,        &
                                   nterms_s     = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%nterms_s,       &
                                   iproc        = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%iproc,          &
                                   itime        = itime_couple,                                                                       &
                                   dof_start    = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%dof_start,      &
                                   recv_comm    = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%recv_comm,      &
                                   recv_domain  = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%recv_domain,    &
                                   recv_element = mesh%domain(idomain_l)%chimera%recv(ChiID)%donor(idepend)%elem_info%recv_element )


                !
                ! Boundary face, linearize wrt bc coupled element
                !
                elseif ( boundary_face ) then

                    group_ID = mesh%domain(idomain_l)%faces(ielement_l,iface)%group_ID
                    patch_ID = mesh%domain(idomain_l)%faces(ielement_l,iface)%patch_ID
                    face_ID  = mesh%domain(idomain_l)%faces(ielement_l,iface)%face_ID

                    call seed%init(idomain_g    = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_g( idepend),  &
                                   idomain_l    = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%idomain_l( idepend),  &
                                   ielement_g   = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_g(idepend),  &
                                   ielement_l   = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%ielement_l(idepend),  &
                                   nfields      = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%nfields(idepend),     &
                                   nterms_s     = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%nterms_s(idepend),    &
                                   iproc        = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%proc(idepend),        &
                                   itime        = itime_couple,                                                                         &
                                   dof_start    = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%dof_start(idepend),   &
                                   recv_comm    = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%recv_comm(idepend),   &
                                   recv_domain  = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%recv_domain(idepend), &
                                   recv_element = mesh%bc_patch_group(group_ID)%patch(patch_ID)%coupling(face_ID)%recv_element(idepend) )


                ! Invalid Case
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
                           nfields      = 0, &
                           nterms_s     = 0, &
                           iproc        = 0, &
                           itime        = 0, &
                           dof_start    = 0, &
                           recv_comm    = 0, &
                           recv_domain  = 0, &
                           recv_element = 0)

        end if


    end function element_compute_seed
    !************************************************************************************************************























end module mod_DNAD_tools
