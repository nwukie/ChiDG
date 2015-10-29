module SCA_boundary_average_advective_flux
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: NFACES,ZERO,ONE,TWO,HALF, &
                                          XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG, &
                                          LOCAL, NEIGHBOR

    use atype_boundary_flux,        only: boundary_flux_t
    use type_mesh,                  only: mesh_t
    use type_solverdata,            only: solverdata_t
    use type_properties,            only: properties_t
    use type_seed,                  only: seed_t
    use type_face_location,         only: face_location_t



    use mod_interpolate,            only: interpolate_face
    use mod_integrate,              only: integrate_volume_flux, integrate_boundary_flux
    use mod_DNAD_tools
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
    !   @param[in]      idom    Domain index
    !   @param[in]      ielem   Element index
    !   @param[in]      iface   Face index
    !   @param[in]      iblk    Block index indicating the linearization direction
    !
    !---------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,idom,ielem,iface,iblk,idonor)
        class(SCA_boundary_average_advective_flux_t),   intent(in)      :: self
        type(mesh_t),                                   intent(in)      :: mesh(:)
        type(solverdata_t),                             intent(inout)   :: sdata
        class(properties_t),                            intent(inout)   :: prop
        integer(ik),                                    intent(in)      :: idom, ielem, iface, iblk
        integer(ik),                                    intent(in)      :: idonor   ! 1 for interior faces, potentially > 1 for Chimera faces

        real(rk)                    :: cx, cy, cz
        integer(ik)                 :: iu, ierr, nnodes
        type(seed_t)                :: seed
        type(face_location_t)       :: face

        type(AD_D), allocatable     :: u_l(:), u_r(:), flux_x(:), flux_y(:), flux_z(:)



        !
        ! Get variable index
        !
        iu        = prop%get_eqn_index('u')


        face%idomain  = idom
        face%ielement = ielem
        face%iface    = iface


        !
        ! Get nGQ nodes
        !
        nnodes    = mesh(idom)%faces(ielem,iface)%gq%nnodes_f



        !
        ! Compute element for linearization
        !
        seed = compute_seed(mesh,idom,ielem,iface,idonor,iblk)



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
        ! Interpolate solution to quadrature nodes
        !
        call interpolate_face(mesh,sdata%q,idom,ielem,iface, iu, u_r, seed, LOCAL)
        call interpolate_face(mesh,sdata%q,idom,ielem,iface, iu, u_l, seed, NEIGHBOR)


        !
        ! Compute boundary average flux
        !
        flux_x = ((cx*u_l + cx*u_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,1)
        flux_y = ((cy*u_l + cy*u_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,2)
        flux_z = ((cz*u_l + cz*u_r)/TWO )  *  mesh(idom)%faces(ielem,iface)%norm(:,3)


        !
        ! Integrate flux
        !
        !call integrate_boundary_flux(mesh(idom)%faces(ielem,iface), sdata, idom, iu, iblk, flux_x, flux_y, flux_z)
        call integrate_boundary_flux(mesh,sdata,face,iu,iblk,idonor,seed,flux_x,flux_y,flux_z)

    end subroutine




end module SCA_boundary_average_advective_flux
