module EULER_volume_advective_flux
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use type_mesh,              only: mesh_t
    use atype_volume_flux,      only: volume_flux_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    
    use mod_interpolate,        only: interpolate_element
    use mod_integrate,          only: integrate_volume_flux
    use mod_DNAD_tools
    use DNAD_D

    use EULER_properties,       only: EULER_properties_t
    implicit none

    private

    
    !> Volume advective flux for Euler equations.
    !!
    !!  @author Nathan A. Wukie
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
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,idom,ielem,iblk)
        class(EULER_volume_advective_flux_t),   intent(in)      :: self
        type(mesh_t),                           intent(in)      :: mesh(:)
        type(solverdata_t),                     intent(inout)   :: sdata
        class(properties_t),                    intent(inout)   :: prop
        integer(ik),                            intent(in)      :: idom, ielem, iblk

        ! Equation indices
        integer(ik)    :: irho
        integer(ik)    :: irhou
        integer(ik)    :: irhov
        integer(ik)    :: irhow
        integer(ik)    :: irhoe


        integer(ik)    :: iseed, i, idonor
        type(seed_t)   :: seed

        type(AD_D), dimension(mesh(idom)%elems(ielem)%gq%vol%nnodes)      ::  &
                    rho, rhou, rhov, rhow, rhoE, p, H,                        &
                    flux_x, flux_y, flux_z


        idonor = 0


        !
        ! Get equation indices
        !
        irho  = prop%get_eqn_index("rho")
        irhou = prop%get_eqn_index("rhou")
        irhov = prop%get_eqn_index("rhov")
        irhow = prop%get_eqn_index("rhow")
        irhoE = prop%get_eqn_index("rhoE")


        !
        ! Get neighbor face and seed element for derivatives
        !
        seed = compute_seed(mesh,idom,ielem,iblk,idonor,iblk)


        !
        ! Interpolate solution to quadrature nodes
        !
        call interpolate_element(mesh,sdata%q,idom,ielem,irho, rho, seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhou,rhou,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhov,rhov,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhow,rhow,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhoE,rhoE,seed)


        !
        ! Compute pressure and total enthalpy
        !
        call prop%fluid%compute_pressure(rho,rhou,rhov,rhow,rhoE,p)

        H = (rhoE + p)/rho

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
        flux_x = (rhou*rhou)/rho  +  p
        flux_y = (rhou*rhov)/rho
        flux_z = (rhou*rhow)/rho

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irhou,iblk,flux_x,flux_y,flux_z)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux_x = (rhov*rhou)/rho
        flux_y = (rhov*rhov)/rho  +  p
        flux_z = (rhov*rhow)/rho

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irhov,iblk,flux_x,flux_y,flux_z)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        flux_x = (rhow*rhou)/rho
        flux_y = (rhow*rhov)/rho
        flux_z = (rhow*rhow)/rho  +  p

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irhow,iblk,flux_x,flux_y,flux_z)

        !============================
        !       ENERGY FLUX
        !============================
        flux_x = rhou*H
        flux_y = rhov*H
        flux_z = rhow*H

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irhoE,iblk,flux_x,flux_y,flux_z)

    end subroutine






end module EULER_volume_advective_flux
