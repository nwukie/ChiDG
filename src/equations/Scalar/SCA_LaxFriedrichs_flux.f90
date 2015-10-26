module SCA_LaxFriedrichs_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ZERO,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use atype_boundary_flux,    only: boundary_flux_t
    use type_mesh,              only: mesh_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    use type_seed,              only: seed_t

    use mod_interpolate,        only: interpolate_face
    use mod_integrate,          only: integrate_boundary_flux
    use mod_DNAD_tools
    use DNAD_D

    use SCA_properties,          only: SCA_properties_t
    implicit none

    private



    !--------------------------------------------------------------------------
    type, extends(boundary_flux_t), public :: SCA_LaxFriedrichs_flux_t


    contains
        procedure   :: compute

    end type SCA_LaxFriedrichs_flux_t
    !---------------------------------------------------------------------------

contains




    !
    !
    !
    !
    !
    !---------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,idom,ielem,iface,iblk,idonor)
        class(SCA_LaxFriedrichs_flux_t),    intent(in)      :: self
        type(mesh_t),                       intent(in)      :: mesh(:)
        type(solverdata_t),                 intent(inout)   :: sdata
        class(properties_t),                intent(inout)   :: prop
        integer(ik),                        intent(in)      :: idom, ielem, iface, iblk
        integer(ik),                        intent(in)      :: idonor   ! 1 for interior faces, potentially > 1 for Chimera faces



        real(rk)                    :: cx, cy, cz
        integer(ik)                 :: iu, ierr, nnodes, ineighbor, i
        integer(ik)                 :: idom_n
        integer(ik)                 :: ielem_n
        integer(ik)                 :: iface_n
        type(seed_t)                :: seed

        type(AD_D), allocatable     :: u_l(:), u_r(:), flux_x(:), flux_y(:), flux_z(:)


        !
        ! Get variable index
        !
        iu        = prop%get_eqn_index('u')

        
        !
        ! Get nGQ nodes
        !
        nnodes    = mesh(idom)%faces(ielem,iface)%gq%nnodes_f


        !
        ! Get neighbor indices
        !
        idom_n    = compute_neighbor_domain( mesh,idom,ielem,iface,idonor)
        ielem_n   = compute_neighbor_element(mesh,idom,ielem,iface,idonor)
        iface_n   = compute_neighbor_face(   mesh,idom,ielem,iface,idonor)


        !
        ! Compute element for linearization
        !
        seed = compute_seed(mesh,idom,ielem,iface,idonor,iblk)
        !idomain_seed = compute_seed_domain( mesh,idom,ielem,iface,iblk,idonor)
        !ielem_seed   = compute_seed_element(mesh,idom,ielem,iface,iblk,idonor)


        !
        ! Get equation set property data
        !
        select type(prop)
            type is (SCA_properties_t)
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
        ! Interpolate solution to quadrature nodes
        !
        call interpolate_face(mesh,sdata%q,idom,   ielem,   iface,  iu,u_r,seed)
        call interpolate_face(mesh,sdata%q,idom_n, ielem_n, iface_n,iu,u_l,seed)


        !
        ! Compute boundary upwind flux
        !
        flux_x = (cx * (u_l - u_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,1)
        flux_y = (cy * (u_l - u_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,2)
        flux_z = (cz * (u_l - u_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,3)


        !
        ! Integrate flux
        !
        call integrate_boundary_flux(mesh(idom)%faces(ielem,iface), sdata, idom, iu, iblk, flux_x, flux_y, flux_z)

    end subroutine









end module SCA_LaxFriedrichs_flux
