module mod_DNAD_tools
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, DIAG
    use type_mesh,              only: mesh_t

    implicit none



contains



    !> Computes the face index of matching face in neighboring element
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  iface   Face index of current element
    !!  @param[out] iface_p Face index of matching face in neighboring element
    !-----------------------------------------------------------------------------
    function compute_neighbor_face(iface) result(iface_p)
        integer(ik), intent(in) :: iface

        integer(ik) :: iface_p

        if (iface == XI_MIN) then
            iface_p = XI_MAX
        else if (iface == XI_MAX) then
            iface_p = XI_MIN
        else if (iface == ETA_MIN) then
            iface_p = ETA_MAX
        else if (iface == ETA_MAX) then
            iface_p = ETA_MIN
        else if (iface == ZETA_MIN) then
            iface_p = ZETA_MAX
        else if (iface == ZETA_MAX) then
            iface_p = ZETA_MIN
        end if

    end function




    !> Compute the index of the element in which we wish to seed derivatives
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  mesh    Mesh definition
    !!  @param[in]  ielem   Index of current element
    !!  @param[in]  iblk    Index of the block linearization we wish to compute
    !!  @param[out] iseed   Index of the element in which we wish to seed derivatives
    !------------------------------------------------------------------------------
    function compute_seed_element(mesh,idom,ielem,iblk) result(iseed)
        type(mesh_t),   intent(in)  :: mesh(:)
        integer(ik),    intent(in)  :: idom
        integer(ik),    intent(in)  :: ielem
        integer(ik),    intent(in)  :: iblk

        integer(ik) :: iseed

        !
        ! Interior Faces
        !
        if (iblk > 0) then

            ! Get element for seeding derivatives
            if ( iblk == DIAG ) then
                iseed = ielem
            else
                iseed = mesh(idom)%faces(ielem,iblk)%ineighbor
            end if


        !
        ! Chimera Faces
        !
        elseif (iblk == 0) then
            call signal(FATAL,"compute_seed_element: derivative seeding has not been implemented for Chimera interfaces yet")
        end if

    end function










end module mod_DNAD_tools
