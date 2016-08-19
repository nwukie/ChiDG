module LA_volume_advective_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF

    use type_mesh,              only: mesh_t
    use type_volume_flux,       only: volume_flux_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    use type_element_info,      only: element_info_t
    use type_function_info,     only: function_info_t

    use mod_interpolate,        only: interpolate
    use mod_integrate,          only: integrate_volume_flux
    use DNAD_D

    use LA_properties,          only: LA_properties_t
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
    type, extends(volume_flux_t), public :: LA_volume_advective_flux_t


    contains
        procedure   :: compute
    end type LA_volume_advective_flux_t
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
        class(LA_volume_advective_flux_t),  intent(in)      :: self
        type(mesh_t),                       intent(in)      :: mesh(:)
        type(solverdata_t),                 intent(inout)   :: sdata
        class(properties_t),                intent(inout)   :: prop
        type(element_info_t),                   intent(in)      :: elem_info
        type(function_info_t),                  intent(in)      :: function_info

        integer(ik)             :: idom, ielem
        type(AD_D), allocatable :: u(:), flux_x(:), flux_y(:), flux_z(:)
        real(rk)                :: cx, cy, cz
        integer(ik)             :: nnodes, ierr
        integer(ik)             :: ivar_u


        idom  = elem_info%idomain_l
        ielem = elem_info%ielement_l


        !
        ! Get variable index from equation set
        !
        ivar_u = prop%get_eqn_index('u')


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
        u = interpolate(mesh,sdata,elem_info,function_info,ivar_u,'value')


        !
        ! Compute volume flux at quadrature nodes
        !
        flux_x = cx  *  u 
        flux_y = cy  *  u
        flux_z = cz  *  u


        !
        ! Integrate volume flux
        !
        call integrate_volume_flux(mesh,sdata,elem_info,function_info,ivar_u,flux_x,flux_y,flux_z)



    end subroutine compute
    !****************************************************************************************************






end module LA_volume_advective_flux
