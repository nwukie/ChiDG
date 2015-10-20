module SCA_volume_advective_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ZERO,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use type_mesh,              only: mesh_t
    use atype_volume_flux,      only: volume_flux_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    use mod_interpolate,        only: interpolate_element
    use mod_integrate,          only: integrate_volume_flux
    use mod_DNAD_tools,         only: compute_neighbor_face, compute_seed_element
    use DNAD_D

    use SCA_properties,          only: SCA_properties_t
    implicit none
    private

    !
    !
    !----------------------------------------------------------------
    type, extends(volume_flux_t), public :: SCA_volume_advective_flux_t


    contains
        procedure   :: compute

    end type SCA_volume_advective_flux_t

contains


    !
    !
    !
    !
    !---------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,idom,ielem,iblk)
        class(SCA_volume_advective_flux_t), intent(in)      :: self
        type(mesh_t),                       intent(in)      :: mesh(:)
        type(solverdata_t),                 intent(inout)   :: sdata
        class(properties_t),                intent(inout)   :: prop
        integer(ik),                        intent(in)      :: idom, ielem, iblk



        type(AD_D), allocatable :: u(:), flux_x(:), flux_y(:), flux_z(:)
        real(rk)                :: cx, cy, cz
        integer(ik)             :: nnodes, ierr, iseed
        integer(ik)             :: ivar_u, i


        associate (elem => mesh(idom)%elems(ielem), q => sdata%q)


            !
            ! Get variable index from equation set
            !
            ivar_u = prop%get_eqn_index('u')


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
            ! Allocate storage for variable values at quadrature points
            !
            nnodes = elem%gq%nnodes_v
            allocate(u(nnodes),         &
                     flux_x(nnodes),    &
                     flux_y(nnodes),    &
                     flux_z(nnodes),    stat = ierr)
            if (ierr /= 0) call AllocationError


            !
            ! Get seed element for derivatives
            !
            iseed   = compute_seed_element(mesh,idom,ielem,iblk)


            !
            ! Interpolate solution to quadrature nodes
            !
            call interpolate_element(mesh,q,idom,ielem,ivar_u,u,iseed)


            !
            ! Compute volume flux at quadrature nodes
            !
            flux_x = cx  *  u 
            flux_y = cy  *  u
            flux_z = cz  *  u


            !
            ! Integrate volume flux
            !
            call integrate_volume_flux(elem,sdata,idom,ivar_u,iblk,flux_x,flux_y,flux_z)

        end associate

    end subroutine






end module SCA_volume_advective_flux
