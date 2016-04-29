module perfect_gas
    use mod_kinds,      only: rk
    use mod_constants,  only: ONE, HALF, ZERO
    use atype_fluid,    only: fluid_t
    use DNAD_D


    !> A fluid type under perfect gas assumptions.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-------------------------------------------------------------------------------------
    type, extends(fluid_t), public :: perfect_gas_t

    contains

        procedure   :: compute_pressure_ad
        procedure   :: compute_pressure_real
        procedure   :: compute_gamma_ad
        procedure   :: compute_gamma_real

    end type perfect_gas_t
    !*************************************************************************************


contains



    !>  Compute pressure of the fluid, under the perfect gas assumption
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]      rho     Fluid density
    !!  @param[in]      rhou    Fluid x-momentum
    !!  @param[in]      rhov    Fluid y-momentum
    !!  @param[in]      rhow    Fluid z-momentum
    !!  @param[in]      rhoE    Fluid Total Energy
    !!  @param[inout]   vals    Fluid pressure
    !!
    !!  TODO: Fix hardcoded gamma
    !!
    !-------------------------------------------------------------------------------------
    subroutine compute_pressure_ad(self,rho,rhou,rhov,rhow,rhoE,vals)
        class(perfect_gas_t),   intent(in)      :: self
        type(AD_D),             intent(in)      :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)
        type(AD_D),             intent(inout)   :: vals(:)

        real(rk)    :: gam

        !& DEBUG, HARDCODED GAMMA
        gam = 1.4_rk

        !vals = (gam-ONE)*(rhoE - HALF*rho*((rhou*rhou)/(rho*rho) + (rhov*rhov)/(rho*rho) + (rhow*rhow)/(rho*rho)))
        vals = (gam-ONE)*(rhoE - HALF*( (rhou*rhou) + (rhov*rhov) + (rhow*rhow) )/rho )

    end subroutine compute_pressure_ad
    !************************************************************************************







    !>  Compute pressure of the fluid, under the perfect gas assumption
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]      rho     Fluid density
    !!  @param[in]      rhou    Fluid x-momentum
    !!  @param[in]      rhov    Fluid y-momentum
    !!  @param[in]      rhow    Fluid z-momentum
    !!  @param[in]      rhoE    Fluid Total Energy
    !!  @param[inout]   vals    Fluid pressure
    !!
    !!  TODO: Fix hardcoded gamma
    !!
    !------------------------------------------------------------------------------------
    subroutine compute_pressure_real(self,rho,rhou,rhov,rhow,rhoE,vals)
        class(perfect_gas_t),   intent(in)      :: self
        real(rk),               intent(in)      :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)
        real(rk),               intent(inout)   :: vals(:)

        real(rk)    :: gam

        !& DEBUG, HARDCODED GAMMA
        gam = 1.4_rk

        vals = (gam-ONE)*(rhoE - HALF*rho*((rhou*rhou)/(rho*rho) + (rhov*rhov)/(rho*rho) + (rhow*rhow)/(rho*rho)))

    end subroutine compute_pressure_real
    !************************************************************************************













    !> Returns a constant gamma value
    !!
    !!
    !!   @author Nathan A. Wukie
    !!   @date   2/1/2016
    !!
    !------------------------------------------------------------------------
    subroutine compute_gamma_ad(self,rho,rhou,rhov,rhow,rhoE,vals)
        class(perfect_gas_t),   intent(in)      :: self
        type(AD_D),             intent(in)      :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)
        type(AD_D),             intent(inout)   :: vals(:)

        !
        ! Make sure vals derivatives are initialized
        !
        vals = ZERO*rho

        !
        ! Set constant value
        !
        vals = 1.4_rk

    end subroutine compute_gamma_ad
    !***********************************************************************




    !> Returns a constant gamma value
    !!
    !!
    !!   @author Nathan A. Wukie
    !!   @date   2/1/2016
    !!
    !------------------------------------------------------------------------
    subroutine compute_gamma_real(self,rho,rhou,rhov,rhow,rhoE,vals)
        class(perfect_gas_t),   intent(in)      :: self
        real(rk),               intent(in)      :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)
        real(rk),               intent(inout)   :: vals(:)

        !
        ! Set constant value
        !
        vals = 1.4_rk

    end subroutine compute_gamma_real
    !***********************************************************************








end module
