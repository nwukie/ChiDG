module SCA_boundary_average_advective_flux
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: NFACES,ZERO,ONE,TWO,HALF, &
                                          XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use atype_boundary_flux,        only: boundary_flux_t
    use type_mesh,                  only: mesh_t
    use atype_solverdata,           only: solverdata_t
    use type_properties,            only: properties_t
    use mod_interpolate,            only: interpolate
    use mod_integrate,              only: integrate_volume_flux, integrate_boundary_flux
    use mod_DNAD_tools,             only: compute_neighbor_face, compute_seed_element
    use DNAD_D

    use SCA_properties,              only: SCA_properties_t
    implicit none

    private

    type, extends(boundary_flux_t), public :: SCA_boundary_average_advective_flux_t


    contains
        procedure   :: compute

    end type SCA_boundary_average_advective_flux_t

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
    subroutine compute(self,mesh,sdata,ielem,iface,iblk,prop)
        class(SCA_boundary_average_advective_flux_t),    intent(in)      :: self
        class(mesh_t),                                  intent(in)      :: mesh
        class(solverdata_t),                            intent(inout)   :: sdata
        integer(ik),                                    intent(in)      :: ielem, iface, iblk
        class(properties_t),                            intent(inout)   :: prop

        real(rk)                    :: cx, cy, cz
        integer(ik)                 :: iu, iseed, ierr, nnodes, ineighbor, iface_p, i
        type(AD_D), allocatable     :: u_l(:), u_r(:), flux_x(:), flux_y(:), flux_z(:)


        iu        = prop%get_eqn_index('u')
        nnodes    = mesh%faces(ielem,iface)%gq%nnodes_f
        ineighbor = mesh%faces(ielem,iface)%ineighbor


        !
        ! Get equation set properties
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
        ! Compute boundary average flux
        !
        flux_x = ((cx*u_l + cx*u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,1)
        flux_y = ((cy*u_l + cy*u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,2)
        flux_z = ((cz*u_l + cz*u_r)/TWO )  *  mesh%faces(ielem,iface)%norm(:,3)


        !
        ! Integrate flux
        !
        call integrate_boundary_flux(mesh%faces(ielem,iface), sdata, iu, iblk, flux_x, flux_y, flux_z)

    end subroutine




end module SCA_boundary_average_advective_flux
