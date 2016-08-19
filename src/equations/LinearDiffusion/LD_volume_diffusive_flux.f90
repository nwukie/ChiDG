module LD_volume_diffusive_flux
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


        type(AD_D), allocatable :: flux_x(:), flux_y(:), flux_z(:), dudx(:), dudy(:), dudz(:)
        real(rk)                :: mu_x, mu_y, mu_z
        integer(ik)             :: iu


        !
        ! Get variable index from equation set
        !
        iu = prop%get_eqn_index('u')


        !
        ! Get equation set properties
        !
        select type(prop)
            type is (LD_properties_t)
                mu_x = prop%mu(1)
                mu_y = prop%mu(2)
                mu_z = prop%mu(3)
        end select


        !
        ! Interpolate solution to quadrature nodes
        !
        dudx = interpolate(mesh,sdata,elem_info,function_info, iu, 'ddx')
        dudy = interpolate(mesh,sdata,elem_info,function_info, iu, 'ddy')
        dudz = interpolate(mesh,sdata,elem_info,function_info, iu, 'ddz')


        !
        ! Compute volume flux at quadrature nodes
        !
        flux_x = mu_x  *  dudx
        flux_y = mu_y  *  dudy
        flux_z = mu_z  *  dudz


        !
        ! Integrate volume flux
        !
        call integrate_volume_flux(mesh,sdata,elem_info,function_info,iu,flux_x,flux_y,flux_z)



    end subroutine compute
    !****************************************************************************************************






end module LD_volume_diffusive_flux
