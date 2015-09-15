module LA_LaxFriedrichs_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ZERO,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use atype_boundary_flux,    only: boundary_flux_t
    use type_mesh,              only: mesh_t
    use atype_solverdata,       only: solverdata_t
    use type_properties,        only: properties_t

    use mod_interpolate,        only: interpolate
    use mod_integrate,          only: integrate_boundary_flux
    use mod_DNAD_tools,         only: compute_neighbor_face, compute_seed_element
    use DNAD_D

    use LA_properties,          only: LA_properties_t
    implicit none

    private



    !--------------------------------------------------------------------------
    type, extends(boundary_flux_t), public :: LA_LaxFriedrichs_flux_t


    contains
        procedure   :: compute

    end type LA_LaxFriedrichs_flux_t
    !---------------------------------------------------------------------------

contains




    !
    !
    !
    !
    !
    !---------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,ielem,iface,iblk,prop)
        class(LA_LaxFriedrichs_flux_t),        intent(in)      :: self
        class(mesh_t),          intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iface, iblk
        class(properties_t),    intent(inout)   :: prop



        real(rk)                 :: cx, cy, cz
        integer(ik)              :: iu, iseed, ierr, nnodes, ineighbor, iface_p, i
        type(AD_D), allocatable  :: u_l(:), u_r(:), flux_x(:), flux_y(:), flux_z(:)


        !
        ! Get integer data
        !
        iu        = prop%get_eqn_index('u')
        nnodes    = mesh%faces(ielem,iface)%gq%nnodes_f
        ineighbor = mesh%faces(ielem,iface)%ineighbor


        !
        ! Get equation set property data
        !
        select type(prop)
            type is (LA_properties_t)
                cx = prop%c(1)
                cy = prop%c(2)
                cz = prop%c(3)
        end select


        !
        ! Allocate arrays for data at quadrature points
        !
        allocate(u_l(nnodes),       &
                 u_r(nnodes),       &
                 flux_x(nnodes),    &
                 flux_y(nnodes),    &
                 flux_z(nnodes), stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Get neighbor face and seed element for derivatives
        !
        iface_p = compute_neighbor_face(iface)
        iseed   = compute_seed_element(mesh,ielem,iblk)


        !
        ! Interpolate solution to quadrature nodes
        !
        call interpolate(mesh%faces,sdata%q,ielem,    iface,  iu,u_r,iseed)
        call interpolate(mesh%faces,sdata%q,ineighbor,iface_p,iu,u_l,iseed)


        !
        ! Compute boundary upwind flux
        !
        flux_x = (cx * (u_l - u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,1)
        flux_y = (cy * (u_l - u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,2)
        flux_z = (cz * (u_l - u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,3)


        !
        ! Integrate flux
        !
        call integrate_boundary_flux(mesh%faces(ielem,iface), sdata, iu, iblk, flux_x, flux_y, flux_z)

    end subroutine









end module LA_LaxFriedrichs_flux
