module LINEULER_volume_advective_source_real
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF,ZERO, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use type_mesh,              only: mesh_t
    use atype_volume_flux,      only: volume_flux_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    
    use mod_interpolate,        only: interpolate_element
    use mod_integrate,          only: integrate_volume_source
    use mod_DNAD_tools
    use DNAD_D

    use LINEULER_properties,    only: LINEULER_properties_t
    implicit none

    private

    type, extends(volume_flux_t), public :: LINEULER_volume_advective_source_real_t


    contains
        procedure  :: compute
    end type LINEULER_volume_advective_source_real_t










contains



    !==========================================================
    !
    !   Volume Flux routine for Euler
    !
    !===========================================================
    subroutine compute(self,mesh,sdata,prop,idom,ielem,iblk)
        class(LINEULER_volume_advective_source_real_t),   intent(in)      :: self
        type(mesh_t),                           intent(in)      :: mesh(:)
        type(solverdata_t),                     intent(inout)   :: sdata
        class(properties_t),                    intent(inout)   :: prop
        integer(ik),                            intent(in)      :: idom, ielem, iblk

        ! Equation indices
        !------------------------------------------------------------
        integer(ik)    :: irho_r, irho_i
        integer(ik)    :: irhou_r, irhou_i
        integer(ik)    :: irhov_r, irhov_i
        integer(ik)    :: irhow_r, irhow_i
        integer(ik)    :: irhoE_r, irhoE_i


        integer(ik)    :: iseed, i, idonor
        type(seed_t)   :: seed

        real(rk)    :: rho_c, rhou_c, rhov_c, rhow_c, rhoE_c
        real(rk)    :: dp_drho, dp_drhou, dp_drhov, dp_drhow, dp_drhoE
        real(rk)    :: ubar, vbar, wbar, Hbar, pbar
        real(rk)    :: gam, omega



        type(AD_D), dimension(mesh(idom)%elems(ielem)%gq%vol%nnodes)      ::  &
                    rho, rhou, rhov, rhow, rhoE, p, H,                        &
                    flux_x, flux_y, flux_z


        idonor = 0


        !-------------------------------------------------------------
        irho_r  = prop%get_eqn_index("rho_r")
        irhou_r = prop%get_eqn_index("rhou_r")
        irhov_r = prop%get_eqn_index("rhov_r")
        irhow_r = prop%get_eqn_index("rhow_r")
        irhoE_r = prop%get_eqn_index("rhoE_r")

        irho_i  = prop%get_eqn_index("rho_i")
        irhou_i = prop%get_eqn_index("rhou_i")
        irhov_i = prop%get_eqn_index("rhov_i")
        irhow_i = prop%get_eqn_index("rhow_i")
        irhoE_i = prop%get_eqn_index("rhoE_i")



        !
        ! Gamma
        !
        gam = 1.4_rk
        omega = 100._rk


        !
        ! Mean flow constants
        !
        rho_c  = 1.2351838930023_rk
        rhou_c = 110.21484155975_rk
        rhov_c = ZERO
        rhow_c = ZERO
        rhoE_c = 267417.20761939_rk


        !
        ! Mean velocities
        !
        ubar = rhou_c / rho_c
        vbar = rhov_c / rho_c
        wbar = rhow_c / rho_c

        !
        ! Mean pressure
        !
        pbar = (gam - ONE) * (rhoE_c - HALF*( (rhou_c*rhou_c) + (rhov_c*rhov_c) + (rhow_c*rhow_c))/rho_c )

        !
        ! Mean enthalpy
        !
        Hbar = (rhoE_c + pbar) / rho_c
        

        !
        ! Pressure jacobians
        !
        dp_drho  = ((gam-ONE)/TWO) * (ubar**TWO + vbar**TWO)
        dp_drhou = -(gam - ONE)*ubar
        dp_drhov = -(gam - ONE)*vbar
        dp_drhow = -(gam - ONE)*wbar
        dp_drhoE = (gam - ONE)









        !
        ! Get neighbor face and seed element for derivatives
        !
        seed = compute_seed(mesh,idom,ielem,iblk,idonor,iblk)




        !
        ! Interpolate solution to quadrature nodes
        !
        call interpolate_element(mesh,sdata%q,idom,ielem,irho_i, rho, seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhou_i,rhou,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhov_i,rhov,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhow_i,rhow,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhoE_i,rhoE,seed)


        !
        ! Compute pressure and total enthalpy
        !
        call prop%fluid%compute_pressure(rho,rhou,rhov,rhow,rhoE,p)

        H = (rhoE + p)/rho

        !===========================
        !        MASS FLUX
        !===========================
        flux_x = -omega * rho

        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,irho_r,iblk,flux_x)


        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux_x = -omega * rhou

        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,irhou_r,iblk,flux_x)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux_x = -omega * rhov

        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,irhov_r,iblk,flux_x)

!        !============================
!        !     Z-MOMENTUM FLUX
!        !============================
!        flux_x = (rhow*rhou)/rho
!        flux_y = (rhow*rhov)/rho
!        flux_z = (rhow*rhow)/rho  +  p
!
!        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irhow,iblk,flux_x,flux_y,flux_z)
!
        !============================
        !       ENERGY FLUX
        !============================
        flux_x = -omega * rhoE

        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,irhoE_r,iblk,flux_x)

    end subroutine






end module LINEULER_volume_advective_source_real
