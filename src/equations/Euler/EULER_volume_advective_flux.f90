module EULER_volume_advective_flux
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF, &
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

    use EULER_properties,       only: EULER_properties_t
    implicit none

    private

    
    !> Volume advective flux for Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(volume_flux_t), public :: EULER_volume_advective_flux_t


    contains

        procedure  :: compute

    end type EULER_volume_advective_flux_t
    !******************************************************************************










contains



    !> Volume flux routine for Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    !subroutine compute(self,mesh,sdata,prop,idom,ielem,iblk)
    subroutine compute(self,mesh,sdata,prop,elem_info,function_info)
        class(EULER_volume_advective_flux_t),   intent(in)      :: self
        type(mesh_t),                           intent(in)      :: mesh(:)
        type(solverdata_t),                     intent(inout)   :: sdata
        class(properties_t),                    intent(inout)   :: prop
        type(element_info_t),                   intent(in)      :: elem_info
        type(function_info_t),                  intent(in)      :: function_info

        ! Equation indices
        integer(ik)    :: irho
        integer(ik)    :: irhou
        integer(ik)    :: irhov
        integer(ik)    :: irhow
        integer(ik)    :: irhoe


        integer(ik)    :: idom, ielem, iblk

        type(AD_D), dimension(mesh(elem_info%idomain_l)%elems(elem_info%ielement_l)%gq%vol%nnodes)      ::  &
                    rho, rhou, rhov, rhow, rhoE, p, H,                        &
                    flux_x, flux_y, flux_z, invrho


        idom  = elem_info%idomain_l
        ielem = elem_info%ielement_l
        iblk  = function_info%iblk


        !
        ! Get equation indices
        !
        irho  = prop%get_eqn_index("rho")
        irhou = prop%get_eqn_index("rhou")
        irhov = prop%get_eqn_index("rhov")
        irhow = prop%get_eqn_index("rhow")
        irhoE = prop%get_eqn_index("rhoE")



        !
        ! Interpolate solution to quadrature nodes
        !
        call interpolate_element(mesh,sdata%q,idom,ielem,irho, rho, function_info%seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhou,rhou,function_info%seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhov,rhov,function_info%seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhow,rhow,function_info%seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhoE,rhoE,function_info%seed)


        invrho = ONE/rho



        !
        ! Compute pressure and total enthalpy
        !
        call prop%fluid%compute_pressure(rho,rhou,rhov,rhow,rhoE,p)

        H = (rhoE + p)*invrho

        !===========================
        !        MASS FLUX
        !===========================
        flux_x = rhou
        flux_y = rhov
        flux_z = rhow

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irho,iblk,flux_x,flux_y,flux_z)


        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux_x = (rhou*rhou)*invrho  +  p
        flux_y = (rhou*rhov)*invrho
        flux_z = (rhou*rhow)*invrho

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irhou,iblk,flux_x,flux_y,flux_z)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux_x = (rhov*rhou)*invrho
        flux_y = (rhov*rhov)*invrho  +  p
        flux_z = (rhov*rhow)*invrho

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irhov,iblk,flux_x,flux_y,flux_z)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        flux_x = (rhow*rhou)*invrho
        flux_y = (rhow*rhov)*invrho
        flux_z = (rhow*rhow)*invrho  +  p

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irhow,iblk,flux_x,flux_y,flux_z)

        !============================
        !       ENERGY FLUX
        !============================
        flux_x = rhou*H
        flux_y = rhov*H
        flux_z = rhow*H

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irhoE,iblk,flux_x,flux_y,flux_z)

    end subroutine compute
    !*********************************************************************************************************






end module EULER_volume_advective_flux
