module DLA_LaxFriedrichs_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ZERO,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use atype_boundary_flux,    only: boundary_flux_t
    use atype_equationset,      only: equationset_t
    use type_mesh,              only: mesh_t
    use atype_solverdata,       only: solverdata_t
    use mod_interpolate,        only: interpolate
    use mod_integrate,          only: integrate_boundary_flux
    use mod_DNAD_tools,         only: compute_neighbor_face, compute_seed_element
    use DNAD_D

    use type_properties,        only: properties_t
    use DLA_properties,         only: DLA_properties_t
    implicit none

    private



    !> This equation set exists really just to test equationsets with 
    !! more than one equation. The idea is just to compute the linear
    !! advecdtion solution twice at the same time. The equations are 
    !! independent of each other. So, we can verify, for example,
    !! the volume flux jacobians for each equation. They should be the
    !! same as for the single LinearAdvection equation set
    !!
    !!
    !-------------------------------------------------------------
    type, extends(boundary_flux_t), public :: DLA_LaxFriedrichs_flux_t

    contains
        procedure   :: compute

    end type DLA_LaxFriedrichs_flux_t

contains

    !==========================================================
    !
    !   Boundary Flux routine for Scalar
    !
    !===========================================================
    subroutine compute(self,mesh,sdata,ielem,iface,iblk,prop)
        class(DLA_LaxFriedrichs_flux_t),    intent(in)      :: self
        class(mesh_t),                      intent(in)      :: mesh
        class(solverdata_t),                intent(inout)   :: sdata
        integer(ik),                        intent(in)      :: ielem, iface, iblk
        class(properties_t),                intent(inout)   :: prop


        real(rk)                 :: cx, cy, cz
        integer(ik)              :: iu_a, iu_b, iseed, ierr, nnodes, ineighbor, iface_p, i
        type(AD_D), allocatable  :: ua_l(:), ua_r(:), ub_l(:), ub_r(:), flux_x(:), flux_y(:), flux_z(:)


        !
        ! Get integer data
        !
        iu_a      = prop%get_eqn_index('u_a')
        iu_b      = prop%get_eqn_index('u_b')
        nnodes    = mesh%faces(ielem,iface)%gq%nnodes_f
        ineighbor = mesh%faces(ielem,iface)%ineighbor

        !
        ! Get equation set property data
        !
        select type(prop)
            type is (DLA_properties_t)
                cx = prop%c(1)
                cy = prop%c(2)
                cz = prop%c(3)
        end select



        !
        ! Allocate arrays for data at quadrature points
        !
        allocate(ua_l(nnodes),       &
                 ua_r(nnodes),       &
                 ub_l(nnodes),       &
                 ub_r(nnodes),       &
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
        call interpolate(mesh%faces,sdata%q,ielem,    iface,  iu_a,ua_r,iseed)
        call interpolate(mesh%faces,sdata%q,ineighbor,iface_p,iu_a,ua_l,iseed)

        call interpolate(mesh%faces,sdata%q,ielem,    iface,  iu_b,ub_r,iseed)
        call interpolate(mesh%faces,sdata%q,ineighbor,iface_p,iu_b,ub_l,iseed)




        !
        ! Compute boundary upwind flux
        !
        flux_x = (cx * (ua_l - ua_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,1)
        flux_y = (cy * (ua_l - ua_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,2)
        flux_z = (cz * (ua_l - ua_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,3)

        call integrate_boundary_flux(mesh%faces(ielem,iface), sdata, iu_a, iblk, flux_x, flux_y, flux_z)



        !> Compute boundary upwind flux
        flux_x = (cx * (ub_l - ub_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,1)
        flux_y = (cy * (ub_l - ub_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,2)
        flux_z = (cz * (ub_l - ub_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,3)

        call integrate_boundary_flux(mesh%faces(ielem,iface), sdata, iu_b, iblk, flux_x, flux_y, flux_z)






    end subroutine








end module DLA_LaxFriedrichs_flux
