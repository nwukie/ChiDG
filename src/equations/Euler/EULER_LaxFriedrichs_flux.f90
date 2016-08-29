module EULER_LaxFriedrichs_flux
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO,HALF,ME,NEIGHBOR

    use type_boundary_flux,     only: boundary_flux_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D

    use EULER_properties,       only: EULER_properties_t
    implicit none

    private


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!
    !----------------------------------------------------------------------------------
    type, extends(boundary_flux_t), public :: EULER_LaxFriedrichs_flux_t

    contains

        procedure  :: compute

    end type EULER_LaxFriedrichs_flux_t
    !**********************************************************************************










contains




    !>  Compute Lax-Friedrichs upwind flux
    !!
    !!  Dissipation = -alpha(u_m - u_p)
    !!
    !!  Alpha is the maximum wave speed
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/16/2016
    !!
    !!------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(EULER_LaxFriedrichs_flux_t),  intent(in)      :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop

        ! Equation indices
        integer(ik)     :: irho
        integer(ik)     :: irhou
        integer(ik)     :: irhov
        integer(ik)     :: irhow
        integer(ik)     :: irhoe


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) :: &
                        rho_m,      rho_p,                                        &
                        rhou_m,     rhou_p,                                       &
                        rhov_m,     rhov_p,                                       &
                        rhow_m,     rhow_p,                                       &
                        rhoe_m,     rhoe_p,                                       &
                        p_m,        p_p,                                          &
                        un_m,       un_p,                                         &
                        a_m,        a_p,                                          &
                        wave_m,     wave_p,                                       &
                        upwind,     wave,                                         &
                        gam_m,      gam_p,                                        &
                        integrand

        real(rk), allocatable, dimension(:)    :: &
                        norm_mag, normx, normy, normz, unormx, unormy, unormz


        irho  = prop%get_eqn_index("rho")
        irhou = prop%get_eqn_index("rhou")
        irhov = prop%get_eqn_index("rhov")
        irhow = prop%get_eqn_index("rhow")
        irhoE = prop%get_eqn_index("rhoE")



        !
        ! Interpolate solution to quadrature nodes
        !
        rho_m  = worker%interpolate(irho,  'value', ME)
        rho_p  = worker%interpolate(irho,  'value', NEIGHBOR)

        rhou_m = worker%interpolate(irhou, 'value', ME)
        rhou_p = worker%interpolate(irhou, 'value', NEIGHBOR)

        rhov_m = worker%interpolate(irhov, 'value', ME)
        rhov_p = worker%interpolate(irhov, 'value', NEIGHBOR)

        rhow_m = worker%interpolate(irhow, 'value', ME)
        rhow_p = worker%interpolate(irhow, 'value', NEIGHBOR)

        rhoE_m = worker%interpolate(irhoE, 'value', ME)
        rhoE_p = worker%interpolate(irhoE, 'value', NEIGHBOR)






        normx  = worker%normal(1)
        normy  = worker%normal(2)
        normz  = worker%normal(3)

        unormx = worker%unit_normal(1)
        unormy = worker%unit_normal(2)
        unormz = worker%unit_normal(3)




        !
        ! Compute pressure and gamma
        !
        call prop%fluid%compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,p_m)
        call prop%fluid%compute_pressure(rho_p,rhou_p,rhov_p,rhow_p,rhoE_p,p_p)
        call prop%fluid%compute_gamma(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,gam_m)
        call prop%fluid%compute_gamma(rho_p,rhou_p,rhov_p,rhow_p,rhoE_p,gam_p)


        !
        ! Compute normal velocities: dot-product vector projection along unit-normal direction
        !
        un_m = unormx*(rhou_m/rho_m) + unormy*(rhov_m/rho_m) + unormz*(rhow_m/rho_m)
        un_p = unormx*(rhou_p/rho_p) + unormy*(rhov_p/rho_p) + unormz*(rhow_p/rho_p)

        
        !
        ! Compute speed of sound
        !
        a_m = sqrt(abs(gam_m * p_m / rho_m))
        a_p = sqrt(abs(gam_p * p_p / rho_p))


        !
        ! Compute wave speeds
        !
        wave_m = abs(un_m) + a_m
        wave_p = abs(un_p) + a_p
        wave   = max(wave_m,wave_p)


        norm_mag = sqrt(normx**TWO + normy**TWO + normz**TWO)

        !================================
        !       MASS FLUX
        !================================
        upwind = -wave*(rho_p - rho_m)

        integrand = HALF * upwind * norm_mag

        call worker%integrate_boundary(irho,integrand)


        !================================
        !       X-MOMENTUM FLUX
        !================================
        upwind = -wave*(rhou_p - rhou_m)

        integrand = HALF * upwind * norm_mag

        call worker%integrate_boundary(irhou,integrand)


        !================================
        !       Y-MOMENTUM FLUX
        !================================
        upwind = -wave*(rhov_p - rhov_m)

        integrand = HALF * upwind * norm_mag

        call worker%integrate_boundary(irhov,integrand)

        !================================
        !       Z-MOMENTUM FLUX
        !================================
        upwind = -wave*(rhow_p - rhow_m)

        integrand = HALF * upwind * norm_mag

        call worker%integrate_boundary(irhow,integrand)

        !================================
        !          ENERGY FLUX
        !================================
        upwind = -wave*(rhoE_p - rhoE_m)

        integrand = HALF * upwind * norm_mag

        call worker%integrate_boundary(irhoE,integrand)


    end subroutine compute
    !*******************************************************************************************













end module EULER_LaxFriedrichs_flux
