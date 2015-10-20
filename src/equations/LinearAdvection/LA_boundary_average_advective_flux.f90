module LA_boundary_average_advective_flux
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: NFACES,ZERO,ONE,TWO,HALF, &
                                          XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use atype_boundary_flux,        only: boundary_flux_t
    use type_mesh,                  only: mesh_t
    use type_solverdata,            only: solverdata_t
    use type_properties,            only: properties_t
    use mod_interpolate,            only: interpolate_face
    use mod_integrate,              only: integrate_boundary_flux
    use mod_DNAD_tools,             only: compute_neighbor_face, compute_seed_element
    use DNAD_D

    use LA_properties,              only: LA_properties_t
    implicit none

    private

    type, extends(boundary_flux_t), public :: LA_boundary_average_advective_flux_t


    contains
        procedure   :: compute

    end type LA_boundary_average_advective_flux_t

contains

    ! Compute the average advective boundary flux for scalar linear advection
    !
    !   @author Nathan A. Wukie
    !
    !   @param[in]      mesh    Mesh data
    !   @param[inout]   sdata   Solver data. Solution, RHS, Linearization etc.
    !   @param[in]      ielem   Element index
    !   @param[in]      iface   Face index
    !   @param[in]      iblk    Block index indicating the linearization direction
    !
    !---------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,idom,ielem,iface,iblk,idonor)
        class(LA_boundary_average_advective_flux_t),    intent(in)      :: self
        type(mesh_t),                                   intent(in)      :: mesh(:)
        type(solverdata_t),                             intent(inout)   :: sdata
        class(properties_t),                            intent(inout)   :: prop
        integer(ik),                                    intent(in)      :: idom, ielem, iface, iblk
        integer(ik),                                    intent(in)      :: idonor

        real(rk)                    :: cx, cy, cz
        integer(ik)                 :: iu, iseed, ierr, nnodes, ineighbor, iface_p, i, idom_n
        type(AD_D), allocatable     :: u_l(:), u_r(:), flux_x(:), flux_y(:), flux_z(:)


        iu        = prop%get_eqn_index('u')
        nnodes    = mesh(idom)%faces(ielem,iface)%gq%nnodes_f
        ineighbor = mesh(idom)%faces(ielem,iface)%ineighbor
        idom_n    = idom


        !
        ! Get equation set properties
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
        iseed   = compute_seed_element(mesh,idom,ielem,iblk)


        !
        ! Interpolate solution to quadrature nodes
        !
        call interpolate_face(mesh,sdata%q,idom,   ielem,    iface,  iu,u_r,iseed)
        call interpolate_face(mesh,sdata%q,idom_n, ineighbor,iface_p,iu,u_l,iseed)


        !
        ! Compute boundary average flux
        !
        flux_x = ((cx*u_l + cx*u_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,1)
        flux_y = ((cy*u_l + cy*u_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,2)
        flux_z = ((cz*u_l + cz*u_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,3)


        !
        ! Integrate flux
        !
        call integrate_boundary_flux(mesh(idom)%faces(ielem,iface), sdata, idom, iu, iblk, flux_x, flux_y, flux_z)

    end subroutine




end module LA_boundary_average_advective_flux
