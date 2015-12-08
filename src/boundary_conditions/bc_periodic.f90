module bc_periodic
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX
    use type_bc,            only: bc_t
    use type_solverdata,    only: solverdata_t
    use type_mesh,          only: mesh_t
    use type_dict,          only: dict_t
    use type_properties,    only: properties_t





    type, extends(bc_t) :: periodic_t

    contains
        procedure init_spec
        procedure compute
    end type periodic_t



contains

    !> Specialized initialization routine called from bc_t%init for matching
    !! periodic boundary conditions. This boundary condition perorms the following:
    !!      - Resets ftype to '0' to indicate an interior face
    !!      - Reset ineighbor to periodic element
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[inout]   mesh    mesh_t containing elements and faces
    !!  @param[inout]   sdata   solverdata_t containing solution(q), rhs, linearization(lin), etc.
    !!  @param[in]      iface   Face index to which the boundary condition is applied
    !!
    !-----------------------------------------------------------------------
    subroutine init_spec(self,mesh,iface,options)
        class(periodic_t),     intent(inout)   :: self
        type(mesh_t),           intent(inout)   :: mesh
        integer(ik),            intent(in)      :: iface
        type(dict_t), optional, intent(in)      :: options

        integer(ik) :: ielem, ielem_p, ixi, ieta, izeta, ixi_p, ieta_p, izeta_p

        call chidg_signal(FATAL,'bc_periodic: broken due to discontinued use of mesh%elem_m matrix-view')

        !
        ! Apply periodic XI
        !
        if (iface == XI_MIN) then
            ixi = 1
            ixi_p = mesh%nelem_xi

            ! Loop through XI-face elements
            do izeta = 1,mesh%nelem_zeta
                do ieta = 1,mesh%nelem_eta
                    ! Get element indices
                    !ielem = mesh%elems_m(ixi,ieta,izeta)%ielem
                    !ielem_p = mesh%elems_m(ixi_p,ieta,izeta)%ielem

                    ! Reset face-type to interior and neighbor to matching periodic element on opposite face
                    mesh%faces(ielem,XI_MIN)%ftype = 0              ! Interior face
                    mesh%faces(ielem,XI_MIN)%ineighbor = ielem_p    ! Set neighbor face to be periodic

                    ! Reset face-type of opposite face and neighbor of opposite face
                    mesh%faces(ielem_p,XI_MAX)%ftype = 0
                    mesh%faces(ielem_p,XI_MAX)%ineighbor = ielem
                end do
            end do
        end if
        

        !
        ! Apply periodic ETA
        !
        if (iface == ETA_MIN) then
            ieta = 1
            ieta_p = mesh%nelem_eta

            ! Loop through ETA-face elements
            do izeta = 1,mesh%nelem_zeta
                do ixi = 1,mesh%nelem_xi
                    ! Get element indices
                    !ielem = mesh%elems_m(ixi,ieta,izeta)%ielem
                    !ielem_p = mesh%elems_m(ixi,ieta_p,izeta)%ielem

                    ! Reset face-type to interior and neighbor to matching periodic element on opposite face
                    mesh%faces(ielem,ETA_MIN)%ftype = 0             ! Interior face
                    mesh%faces(ielem,ETA_MIN)%ineighbor = ielem_p   ! Set neighbor face to be periodic

                    ! Reset face-type of opposite face and neighbor of opposite face
                    mesh%faces(ielem_p,ETA_MAX)%ftype = 0
                    mesh%faces(ielem_p,ETA_MAX)%ineighbor = ielem
                end do
            end do
        end if


        !
        ! Apply periodic ZETA
        !
        if (iface == ZETA_MIN) then
            izeta = 1
            izeta_p = mesh%nelem_zeta

            ! Loop through ZETA-face elements
            do ieta = 1,mesh%nelem_eta
                do ixi = 1,mesh%nelem_xi
                    ! Get element indices
                    !ielem = mesh%elems_m(ixi,ieta,izeta)%ielem
                    !ielem_p = mesh%elems_m(ixi,ieta,izeta_p)%ielem

                    ! Reset face-type to interior and neighbor to matching periodic element on opposite face
                    mesh%faces(ielem,ZETA_MIN)%ftype = 0            ! Interior face
                    mesh%faces(ielem,ZETA_MIN)%ineighbor = ielem_p  ! Set neighbor face to be periodic

                    ! Reset face-type of opposite face and neighbor of opposite face
                    mesh%faces(ielem_p,ZETA_MAX)%ftype = 0
                    mesh%faces(ielem_p,ZETA_MAX)%ineighbor = ielem
                end do
            end do
        end if



    end subroutine init_spec







    !> Boundary condition compute routine called by spatial scheme
    !!      - Matching periodic boundary condition, so the interior scheme 
    !!        is just adjusted to connect elements and no extra calculation
    !!        routine is needed here. Hence the empty routine below.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!
    !------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,idom,ielem,iface,iblk)
        class(periodic_t),      intent(inout)   :: self
        type(mesh_t),           intent(in)      :: mesh(:)
        type(solverdata_t),     intent(inout)   :: sdata
        class(properties_t),    intent(inout)   :: prop
        integer(ik),            intent(in)      :: idom
        integer(ik),            intent(in)      :: ielem
        integer(ik),            intent(in)      :: iface
        integer(ik),            intent(in)      :: iblk

        ! DO NOTHING IN PERIODIC BOUNDARY CONDITION

    end subroutine compute













end module
