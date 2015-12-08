module mod_DNAD_tools
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, DIAG, CHIMERA, INTERIOR, BOUNDARY
    use type_mesh,              only: mesh_t
    use type_seed,              only: seed_t
    implicit none



contains










    !> Computes the domain index of the neighbor domain
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  mesh    Array of mesh_t instances
    !!  @param[in]  idom    Domain index of the current element
    !!  @param[in]  ielem   Element index for mesh(idom)%elems(ielem)
    !!  @param[in]  iface   Face index for mesh(idom)%faces(ielem,iface)
    !!  @param[in]  idonor  Index of potential donor elements. For example, Chimera faces could have more than one donor element
    !!
    !-----------------------------------------------------------------------------------------------
    function compute_neighbor_domain(mesh,idom,ielem,iface,idonor) result(idom_n)
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
            idom_n = mesh(idom)%chimera%recv%data(ChiID)%donor_domain%at(idonor)

        else
            idom_n = idom
        end if

    end function









    !> Computes the element index of the neighbor element
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  mesh    Array of mesh_t instances
    !!  @param[in]  idom    Domain index of the current element
    !!  @param[in]  ielem   Element index for mesh(idom)%elems(ielem)
    !!  @param[in]  iface   Face index for mesh(idom)%faces(ielem,iface)
    !!  @param[in]  idonor  Index of potential donor elements. For example, Chimera faces could have more than one donor element
    !!
    !------------------------------------------------------------------------------------------------
    function compute_neighbor_element(mesh,idom,ielem,iface,idonor) result(ielem_n)
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
            ielem_n = mesh(idom)%chimera%recv%data(ChiID)%donor_element%at(idonor)

        else
            ielem_n = mesh(idom)%faces(ielem,iface)%ineighbor
        end if


    end function






    !> Computes the face index of matching face in neighboring element
    !!
    !!  @author Nathan A. Wukie
    !!
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

    end function












    !> Compute the domain and element that is being linearized. These are found by checking
    !! the current element/face neighbors, or by setting the current element itself. Which
    !! element gets linearized depends on iblk. iblk specifies the direction of the linearization.
    !! For example, if iblk == XI_MIN, the seed element is that which is the neighbor to the XI_MIN
    !! face of the current element. The other behavior is if the linearization of the current
    !! element is desired. That is handled with iblk == DIAG.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  mesh    Array of mesh_t instances
    !!  @param[in]  idom    Domain index of the current element
    !!  @param[in]  ielem   Element index for mesh(idom)%elems(ielem)
    !!  @param[in]  iface   Face index for mesh(idom)%faces(ielem,iface)
    !!  @param[in]  idonor  Index of potential donor elements. For example, Chimera faces could have more than one donor element
    !!  @param[in]  iblk    Linearization index
    !!
    !-------------------------------------------------------------------------------------------------------------------------
    function compute_seed(mesh,idom,ielem,iface,idonor,iblk) result(seed)
        type(mesh_t),   intent(in)  :: mesh(:)
        integer(ik),    intent(in)  :: idom
        integer(ik),    intent(in)  :: ielem
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

            seed%idom    = idom
            seed%ielem   = ielem



        !
        ! Check for linearization of neighbor element (neighbor of face(ielem,iface) )
        !
        elseif ( linearize_neighbor ) then



            !
            ! Check if linearization direction (iface) is a Interior or Chimera face
            !
            chimera_face  = ( mesh(idom)%faces(ielem,iface)%ftype == CHIMERA  )
            interior_face = ( mesh(idom)%faces(ielem,iface)%ftype == INTERIOR )
            boundary_face = ( mesh(idom)%faces(ielem,iface)%ftype == BOUNDARY )


            !
            ! Linearize wrt Standard Interior Neighbor
            !
            if ( interior_face ) then

                seed%idom  = idom
                seed%ielem = mesh(idom)%faces(ielem,iface)%ineighbor


            !
            ! Linearize wrt Chimera Interior Neighbor
            !
            elseif ( chimera_face ) then
                ChiID = mesh(idom)%faces(ielem,iface)%ChiID

                seed%idom  = mesh(idom)%chimera%recv%data(ChiID)%donor_domain%at(idonor)
                seed%ielem = mesh(idom)%chimera%recv%data(ChiID)%donor_element%at(idonor)


            !
            ! Boudnary face, linearize wrt Current element
            !
            elseif ( boundary_face ) then
                seed%idom  = idom
                seed%ielem = ielem

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















end module mod_DNAD_tools
