module mod_DNAD_tools
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, DIAG, CHIMERA, INTERIOR, BOUNDARY
    use mod_chidg_mpi,          only: IRANK
    
    use type_mesh,              only: mesh_t
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
        type(mesh_t),   intent(in)  :: mesh(:)
        integer(ik),    intent(in)  :: idom
        integer(ik),    intent(in)  :: ielem
        integer(ik),    intent(in)  :: iface
        integer(ik),    intent(in)  :: idonor

        integer(ik) :: idom_n, ChiID
        logical     :: chimera_face = .false.

        chimera_face = (mesh(idom)%faces(ielem,iface)%ftype == CHIMERA)


        if ( chimera_face ) then

            ChiID  = mesh(idom)%faces(ielem,iface)%ChiID
            idom_n = mesh(idom)%chimera%recv%data(ChiID)%donor_domain_l%at(idonor)

        else
            !idom_n = idom
            idom_n = mesh(idom)%faces(ielem,iface)%ineighbor_domain_l
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
        type(mesh_t),   intent(in)  :: mesh(:)
        integer(ik),    intent(in)  :: idom
        integer(ik),    intent(in)  :: ielem
        integer(ik),    intent(in)  :: iface
        integer(ik),    intent(in)  :: idonor

        integer(ik) :: idom_n, ChiID
        logical     :: chimera_face = .false.

        chimera_face = (mesh(idom)%faces(ielem,iface)%ftype == CHIMERA)


        if ( chimera_face ) then

            ChiID  = mesh(idom)%faces(ielem,iface)%ChiID
            idom_n = mesh(idom)%chimera%recv%data(ChiID)%donor_domain_g%at(idonor)

        else
            !idom_n = idom
            idom_n = mesh(idom)%faces(ielem,iface)%ineighbor_domain_g
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
        type(mesh_t),   intent(in)  :: mesh(:)
        integer(ik),    intent(in)  :: idom
        integer(ik),    intent(in)  :: ielem
        integer(ik),    intent(in)  :: iface
        integer(ik),    intent(in)  :: idonor

        integer(ik) :: ielem_n, ChiID
        logical     :: chimera_face = .false.

        chimera_face = (mesh(idom)%faces(ielem,iface)%ftype == CHIMERA)
        
        if ( chimera_face ) then

            ChiID   = mesh(idom)%faces(ielem,iface)%ChiID
            ielem_n = mesh(idom)%chimera%recv%data(ChiID)%donor_element_l%at(idonor)

        else
            ielem_n = mesh(idom)%faces(ielem,iface)%get_neighbor_element_l()
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
        type(mesh_t),   intent(in)  :: mesh(:)
        integer(ik),    intent(in)  :: idom
        integer(ik),    intent(in)  :: ielem
        integer(ik),    intent(in)  :: iface
        integer(ik),    intent(in)  :: idonor

        integer(ik) :: ielem_n, ChiID
        logical     :: chimera_face = .false.

        chimera_face = (mesh(idom)%faces(ielem,iface)%ftype == CHIMERA)
        
        if ( chimera_face ) then

            ChiID   = mesh(idom)%faces(ielem,iface)%ChiID
            ielem_n = mesh(idom)%chimera%recv%data(ChiID)%donor_element_l%at(idonor)

        else
            ielem_n = mesh(idom)%faces(ielem,iface)%ineighbor_element_g
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
        type(mesh_t),   intent(in)  :: mesh(:)
        integer(ik),    intent(in)  :: idom
        integer(ik),    intent(in)  :: ielem
        integer(ik),    intent(in)  :: iface
        integer(ik),    intent(in)  :: idonor

        integer(ik) :: iface_n
        logical     :: chimera_face

        chimera_face = ( mesh(idom)%faces(ielem,iface)%ftype == CHIMERA )

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
    !! element gets linearized depends on iblk. iblk specifies the direction of the linearization.
    !! For example, if iblk == XI_MIN, the seed element is that which is the neighbor to the XI_MIN
    !! face of the current element. The other behavior is if the linearization of the current
    !! element is desired. That is handled with iblk == DIAG.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  mesh    Array of mesh_t instances
    !!  @param[in]  idom    Domain index of the current element
    !!  @param[in]  ielem   Element index for mesh(idom)%elems(ielem)
    !!  @param[in]  iface   Face index for mesh(idom)%faces(ielem,iface)
    !!  @param[in]  idonor  Index of potential donor elements. For example, Chimera faces could have more than one donor element
    !!  @param[in]  iblk    Linearization index
    !!
    !-------------------------------------------------------------------------------------------------------------------------
    function compute_seed(mesh,idomain_l,ielement_l,iface,idonor,iblk) result(seed)
        type(mesh_t),   intent(in)  :: mesh(:)
        integer(ik),    intent(in)  :: idomain_l
        integer(ik),    intent(in)  :: ielement_l
        integer(ik),    intent(in)  :: iface
        integer(ik),    intent(in)  :: idonor
        integer(ik),    intent(in)  :: iblk


        type(seed_t)    :: seed
        integer(ik)     :: ChiID
        logical         :: linearize_me, linearize_neighbor
        logical         :: chimera_face, interior_face, boundary_face



        linearize_me        = ( iblk == DIAG )
        linearize_neighbor  = ( iblk == XI_MIN   .or. iblk == XI_MAX   .or. &
                                iblk == ETA_MIN  .or. iblk == ETA_MAX  .or. &
                                iblk == ZETA_MIN .or. iblk == ZETA_MAX )



        !
        ! Check for linearization of current element (ielem)
        !
        if ( linearize_me ) then

            seed%idomain_g  = mesh(idomain_l)%elems(ielement_l)%idomain_g
            seed%idomain_l  = mesh(idomain_l)%elems(ielement_l)%idomain_l
            seed%ielement_g = mesh(idomain_l)%elems(ielement_l)%ielement_g
            seed%ielement_l = mesh(idomain_l)%elems(ielement_l)%ielement_l
            seed%iproc      = IRANK


        !
        ! Check for linearization of neighbor element (neighbor of face(ielem,iface) )
        !
        elseif ( linearize_neighbor ) then



            !
            ! Check if linearization direction (iface) is a Interior or Chimera face
            !
            chimera_face  = ( mesh(idomain_l)%faces(ielement_l,iface)%ftype == CHIMERA  )
            interior_face = ( mesh(idomain_l)%faces(ielement_l,iface)%ftype == INTERIOR )
            boundary_face = ( mesh(idomain_l)%faces(ielement_l,iface)%ftype == BOUNDARY )


            !
            ! Linearize wrt Standard Interior Neighbor
            !
            if ( interior_face ) then

                seed%idomain_g    = mesh(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_g
                seed%idomain_l    = mesh(idomain_l)%faces(ielement_l,iface)%ineighbor_domain_l
                seed%ielement_g   = mesh(idomain_l)%faces(ielement_l,iface)%ineighbor_element_g
                seed%ielement_l   = mesh(idomain_l)%faces(ielement_l,iface)%ineighbor_element_l
                seed%iproc        = mesh(idomain_l)%faces(ielement_l,iface)%ineighbor_proc
                seed%recv_comm    = mesh(idomain_l)%faces(ielement_l,iface)%recv_comm
                seed%recv_domain  = mesh(idomain_l)%faces(ielement_l,iface)%recv_domain
                seed%recv_element = mesh(idomain_l)%faces(ielement_l,iface)%recv_element


            !
            ! Linearize wrt Chimera Interior Neighbor
            !
            elseif ( chimera_face ) then
                ChiID = mesh(idomain_l)%faces(ielement_l,iface)%ChiID

                seed%idomain_g    = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_domain_g%at(idonor)
                seed%idomain_l    = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_domain_l%at(idonor)
                seed%ielement_g   = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_element_g%at(idonor)
                seed%ielement_l   = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_element_l%at(idonor)
                seed%iproc        = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_proc%at(idonor)
                seed%recv_comm    = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_recv_comm%at(idonor)
                seed%recv_domain  = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_recv_domain%at(idonor)
                seed%recv_element = mesh(idomain_l)%chimera%recv%data(ChiID)%donor_recv_element%at(idonor)


            !
            ! Boudnary face, linearize wrt Current element
            !
            elseif ( boundary_face ) then
                seed%idomain_g  = mesh(idomain_l)%elems(ielement_l)%idomain_g
                seed%idomain_l  = mesh(idomain_l)%elems(ielement_l)%idomain_l
                seed%ielement_g = mesh(idomain_l)%elems(ielement_l)%ielement_g
                seed%ielement_l = mesh(idomain_l)%elems(ielement_l)%ielement_l
                seed%iproc      = IRANK


            !
            ! Invalid Case
            !
            else
                call chidg_signal(FATAL,"compute_seed_domain: invalid face type - face(ielem,iface)%ftype")
            end if


        !
        ! Invalid value for iface
        !
        else
            call chidg_signal(FATAL,"compute_seed_domain: invalid value for iface")

        end if





    end function compute_seed
    !************************************************************************************************************















end module mod_DNAD_tools
