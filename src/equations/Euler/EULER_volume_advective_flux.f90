module EULER_volume_advective_flux
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX

    use atype_volume_flux,      only: volume_flux_t
    use atype_equationset,      only: equationset_t
    use type_mesh,              only: mesh_t
    use atype_solverdata,       only: solverdata_t
    use mod_interpolate,        only: interpolate
    use mod_integrate,          only: integrate_volume_flux
    use mod_DNAD_tools,         only: compute_neighbor_face, compute_seed_element
    use DNAD_D

    use type_properties,        only: properties_t
    use EULER_properties,       only: EULER_properties_t
    implicit none

    private

    type, extends(volume_flux_t), public :: EULER_volume_advective_flux_t


    contains
        procedure  :: compute
    end type EULER_volume_advective_flux_t










contains



    !==========================================================
    !
    !   Volume Flux routine for Euler
    !
    !===========================================================
    subroutine compute(self,mesh,sdata,ielem,iblk,prop)
        class(EULER_volume_advective_flux_t),   intent(in)      :: self
        class(mesh_t),                          intent(in)      :: mesh
        class(solverdata_t),                    intent(inout)   :: sdata
        integer(ik),                            intent(in)      :: ielem, iblk
        class(properties_t),                    intent(inout)   :: prop

        ! Equation indices
        !------------------------------------------------------------
        integer(ik)    :: irho
        integer(ik)    :: irhou
        integer(ik)    :: irhov
        integer(ik)    :: irhow
        integer(ik)    :: irhoe


        integer(ik)    :: iseed, i

        type(AD_D), dimension(mesh%elems(ielem)%gq%vol%nnodes)      ::  &
                    rho, rhou, rhov, rhow, rhoE, p, H,                  &
                    flux_x, flux_y, flux_z

        !-------------------------------------------------------------
        irho  = prop%get_eqn_index("rho")
        irhou = prop%get_eqn_index("rhou")
        irhov = prop%get_eqn_index("rhov")
        irhow = prop%get_eqn_index("rhow")
        irhoE = prop%get_eqn_index("rhoE")

        !> Get neighbor face and seed element for derivatives
        iseed   = compute_seed_element(mesh,ielem,iblk)


        !> Interpolate solution to quadrature nodes
        call interpolate(mesh%elems,sdata%q,ielem,irho, rho, iseed)
        call interpolate(mesh%elems,sdata%q,ielem,irhou,rhou,iseed)
        call interpolate(mesh%elems,sdata%q,ielem,irhov,rhov,iseed)
        call interpolate(mesh%elems,sdata%q,ielem,irhow,rhow,iseed)
        call interpolate(mesh%elems,sdata%q,ielem,irhoE,rhoE,iseed)


        ! Compute pressure and total enthalpy
        call prop%fluid%compute_pressure(rho,rhou,rhov,rhow,rhoE,p)

        H = (rhoE + p)/rho

        !===========================
        !        MASS FLUX
        !===========================
        flux_x = rhou
        flux_y = rhov
        flux_z = rhow

        call integrate_volume_flux(mesh%elems(ielem),sdata,irho,iblk,flux_x,flux_y,flux_z)


        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux_x = (rhou*rhou)/rho  +  p
        flux_y = (rhou*rhov)/rho
        flux_z = (rhou*rhow)/rho

        call integrate_volume_flux(mesh%elems(ielem),sdata,irhou,iblk,flux_x,flux_y,flux_z)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux_x = (rhov*rhou)/rho
        flux_y = (rhov*rhov)/rho  +  p
        flux_z = (rhov*rhow)/rho

        call integrate_volume_flux(mesh%elems(ielem),sdata,irhov,iblk,flux_x,flux_y,flux_z)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        flux_x = (rhow*rhou)/rho
        flux_y = (rhow*rhov)/rho
        flux_z = (rhow*rhow)/rho  +  p

        call integrate_volume_flux(mesh%elems(ielem),sdata,irhow,iblk,flux_x,flux_y,flux_z)

        !============================
        !       ENERGY FLUX
        !============================
        flux_x = rhou*H
        flux_y = rhov*H
        flux_z = rhow*H

        call integrate_volume_flux(mesh%elems(ielem),sdata,irhoE,iblk,flux_x,flux_y,flux_z)

    end subroutine






end module EULER_volume_advective_flux
