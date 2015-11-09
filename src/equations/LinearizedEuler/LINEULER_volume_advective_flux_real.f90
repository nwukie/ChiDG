module LINEULER_volume_advective_flux_real
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF,ZERO, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use type_mesh,              only: mesh_t
    use atype_volume_flux,      only: volume_flux_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    
    use mod_interpolate,        only: interpolate_element
    use mod_integrate,          only: integrate_volume_flux
    use mod_DNAD_tools
    use DNAD_D

    use LINEULER_properties,    only: LINEULER_properties_t
    use mod_linearized_euler
    implicit none

    private

    type, extends(volume_flux_t), public :: LINEULER_volume_advective_flux_real_t


    contains
        procedure  :: compute
    end type LINEULER_volume_advective_flux_real_t










contains



    !==========================================================
    !
    !   Volume Flux routine for Euler
    !
    !===========================================================
    subroutine compute(self,mesh,sdata,prop,idom,ielem,iblk)
        class(LINEULER_volume_advective_flux_real_t),   intent(in)      :: self
        type(mesh_t),                           intent(in)      :: mesh(:)
        type(solverdata_t),                     intent(inout)   :: sdata
        class(properties_t),                    intent(inout)   :: prop
        integer(ik),                            intent(in)      :: idom, ielem, iblk

        ! Equation indices
        !------------------------------------------------------------
        integer(ik)    :: irho
        integer(ik)    :: irhou
        integer(ik)    :: irhov
        integer(ik)    :: irhow
        integer(ik)    :: irhoe


        integer(ik)    :: iseed, i, idonor
        type(seed_t)   :: seed




        type(AD_D), dimension(mesh(idom)%elems(ielem)%gq%vol%nnodes)      ::  &
                    rho, rhou, rhov, rhow, rhoE,                        &
                    flux_x, flux_y, flux_z


        idonor = 0


        !-------------------------------------------------------------
        irho  = prop%get_eqn_index("rho_r")
        irhou = prop%get_eqn_index("rhou_r")
        irhov = prop%get_eqn_index("rhov_r")
        irhow = prop%get_eqn_index("rhow_r")
        irhoE = prop%get_eqn_index("rhoE_r")





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



        !===========================
        !        MASS FLUX
        !===========================
        flux_x = rho_x_rho  * rho  + &
                 rho_x_rhou * rhou + &
                 rho_x_rhov * rhov + &
                 rho_x_rhoE * rhoE
        flux_y = rho_y_rho  * rho  + &
                 rho_y_rhou * rhou + &
                 rho_y_rhov * rhov + &
                 rho_y_rhoE * rhoE
        flux_z = rhow

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irho,iblk,flux_x,flux_y,flux_z)


        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux_x = rhou_x_rho  * rho  + &
                 rhou_x_rhou * rhou + &
                 rhou_x_rhov * rhov + &
                 rhou_x_rhoE * rhoE
        flux_y = rhou_y_rho  * rho  + &
                 rhou_y_rhou * rhou + &
                 rhou_y_rhov * rhov + &
                 rhou_y_rhoE * rhoE
        flux_z = ZERO

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irhou,iblk,flux_x,flux_y,flux_z)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux_x = rhov_x_rho  * rho  + &
                 rhov_x_rhou * rhou + &
                 rhov_x_rhov * rhov + &
                 rhov_x_rhoE * rhoE
        flux_y = rhov_y_rho  * rho  + &
                 rhov_y_rhou * rhou + &
                 rhov_y_rhov * rhov + &
                 rhov_y_rhoE * rhoE
        flux_z = ZERO

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irhov,iblk,flux_x,flux_y,flux_z)

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
        flux_x = rhoE_x_rho  * rho  + &
                 rhoE_x_rhou * rhou + &
                 rhoE_x_rhov * rhov + &
                 rhoE_x_rhoE * rhoE
        flux_y = rhoE_y_rho  * rho  + &
                 rhoE_y_rhou * rhou + &
                 rhoE_y_rhov * rhov + &
                 rhoE_y_rhoE * rhoE
        flux_z = ZERO

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irhoE,iblk,flux_x,flux_y,flux_z)

    end subroutine






end module LINEULER_volume_advective_flux_real
