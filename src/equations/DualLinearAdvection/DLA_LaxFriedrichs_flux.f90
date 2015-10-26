module DLA_LaxFriedrichs_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ZERO,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use atype_boundary_flux,    only: boundary_flux_t
    use type_mesh,              only: mesh_t
    use type_solverdata,        only: solverdata_t
    use type_seed,              only: seed_t


    use mod_interpolate,        only: interpolate_face
    use mod_integrate,          only: integrate_boundary_flux
    use mod_DNAD_tools,         only: compute_neighbor_face, compute_seed
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
    subroutine compute(self,mesh,sdata,prop,idom,ielem,iface,iblk,idonor)
        class(DLA_LaxFriedrichs_flux_t),    intent(in)      :: self
        type(mesh_t),                       intent(in)      :: mesh(:)
        type(solverdata_t),                 intent(inout)   :: sdata
        class(properties_t),                intent(inout)   :: prop
        integer(ik),                        intent(in)      :: idom, ielem, iface, iblk
        integer(ik),                        intent(in)      :: idonor


        real(rk)                 :: cx, cy, cz
        integer(ik)              :: iu_a, iu_b, ierr, nnodes, ineighbor, iface_p, i, idom_n
        type(seed_t)             :: seed
        type(AD_D), allocatable  :: ua_l(:), ua_r(:), ub_l(:), ub_r(:), flux_x(:), flux_y(:), flux_z(:)


        !
        ! Get integer data
        !
        iu_a      = prop%get_eqn_index('u_a')
        iu_b      = prop%get_eqn_index('u_b')
        nnodes    = mesh(idom)%faces(ielem,iface)%gq%nnodes_f
        ineighbor = mesh(idom)%faces(ielem,iface)%ineighbor
        idom_n    = idom

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
        iface_p = compute_neighbor_face(mesh,idom,ielem,iface,idonor)


        !
        ! Compute element for linearization
        !
        seed = compute_seed(mesh,idom,ielem,iface,idonor,iblk)
        !iseed   = compute_seed_element(mesh,idom,ielem,iface,iblk,idonor)


        !
        ! Interpolate solution to quadrature nodes
        !
        call interpolate_face(mesh,sdata%q,idom,   ielem,    iface,  iu_a,ua_r,seed)
        call interpolate_face(mesh,sdata%q,idom_n, ineighbor,iface_p,iu_a,ua_l,seed)

        call interpolate_face(mesh,sdata%q,idom,   ielem,    iface,  iu_b,ub_r,seed)
        call interpolate_face(mesh,sdata%q,idom_n, ineighbor,iface_p,iu_b,ub_l,seed)




        !
        ! Compute boundary upwind flux
        !
        flux_x = (cx * (ua_l - ua_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,1)
        flux_y = (cy * (ua_l - ua_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,2)
        flux_z = (cz * (ua_l - ua_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,3)

        call integrate_boundary_flux(mesh(idom)%faces(ielem,iface), sdata, idom, iu_a, iblk, flux_x, flux_y, flux_z)



        !> Compute boundary upwind flux
        flux_x = (cx * (ub_l - ub_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,1)
        flux_y = (cy * (ub_l - ub_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,2)
        flux_z = (cz * (ub_l - ub_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,3)

        call integrate_boundary_flux(mesh(idom)%faces(ielem,iface), sdata, idom, iu_b, iblk, flux_x, flux_y, flux_z)






    end subroutine








end module DLA_LaxFriedrichs_flux
