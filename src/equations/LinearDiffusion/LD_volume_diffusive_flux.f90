module LD_volume_diffusive_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ZERO,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use type_mesh,              only: mesh_t
    use type_volume_flux,       only: volume_flux_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    use type_element_info,      only: element_info_t
    use type_function_info,     only: function_info_t

    use mod_interpolate,        only: interpolate_element
    use mod_integrate,          only: integrate_volume_flux
    use DNAD_D

    use LD_properties,          only: LD_properties_t
    implicit none
    private

    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------
    type, extends(volume_flux_t), public :: LD_volume_diffusive_flux_t


    contains
        procedure   :: compute
    end type LD_volume_diffusive_flux_t
    !*************************************************************************

contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,elem_info,function_info)
        class(LD_volume_diffusive_flux_t),  intent(in)      :: self
        type(mesh_t),                       intent(in)      :: mesh(:)
        type(solverdata_t),                 intent(inout)   :: sdata
        class(properties_t),                intent(inout)   :: prop
        type(element_info_t),               intent(in)      :: elem_info
        type(function_info_t),              intent(in)      :: function_info



        integer(ik)             :: idom, ielem, iblk
        type(AD_D), allocatable :: u(:), flux_x(:), flux_y(:), flux_z(:)
        real(rk)                :: cx, cy, cz
        integer(ik)             :: nnodes, ierr
        integer(ik)             :: ivar_u


        idom  = elem_info%idomain_l
        ielem = elem_info%ielement_l
        iblk  = function_info%iblk


        !
        ! Get variable index from equation set
        !
        ivar_u = prop%get_eqn_index('u')


        !
        ! Get equation set properties
        !
        select type(prop)
            type is (LD_properties_t)
                cx = prop%c(1)
                cy = prop%c(2)
                cz = prop%c(3)
        end select


        !
        ! Allocate storage for variable values at quadrature points
        !
        nnodes = mesh(idom)%elems(ielem)%gq%nnodes_v

        
        allocate(u(nnodes),         &
                 flux_x(nnodes),    &
                 flux_y(nnodes),    &
                 flux_z(nnodes),    stat = ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Interpolate solution to quadrature nodes
        !
        call interpolate_element(mesh,sdata%q,idom,ielem,ivar_u,u,function_info%seed)


        !
        ! Compute volume flux at quadrature nodes
        !
        flux_x = cx  *  u 
        flux_y = cy  *  u
        flux_z = cz  *  u


        !
        ! Integrate volume flux
        !
        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,ivar_u,iblk,flux_x,flux_y,flux_z)



    end subroutine compute
    !****************************************************************************************************






end module LD_volume_diffusive_flux
