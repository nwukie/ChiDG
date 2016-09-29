module perfect_gas
    use mod_kinds,      only: rk
    use mod_constants,  only: ZERO, HALF, ONE, TWO, THREE
    use type_fluid,     only: fluid_t
    use DNAD_D
    implicit none


    !> A fluid type under perfect gas assumptions.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-------------------------------------------------------------------------------------
    type, extends(fluid_t), public :: perfect_gas_t

        real(rk) :: R = 287.15_rk

    contains

        procedure   :: compute_pressure_ad
        procedure   :: compute_pressure_real
        procedure   :: compute_gamma_ad
        procedure   :: compute_gamma_real
        procedure   :: compute_temperature_ad
        procedure   :: compute_viscosity_dynamic_ad
        procedure   :: compute_viscosity_second_ad

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
    function compute_pressure_ad(self,rho,rhou,rhov,rhow,rhoE) result(vals)
        class(perfect_gas_t),           intent(in)      :: self
        type(AD_D),                     intent(in)      :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)

        type(AD_D),     allocatable :: vals(:)

        real(rk) :: gam

        !& DEBUG, HARDCODED GAMMA
        gam = 1.4_rk

        vals = (gam-ONE)*(rhoE - HALF*( (rhou*rhou) + (rhov*rhov) + (rhow*rhow) )/rho )

    end function compute_pressure_ad
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
    function compute_pressure_real(self,rho,rhou,rhov,rhow,rhoE) result(vals)
        class(perfect_gas_t),       intent(in)      :: self
        real(rk),                   intent(in)      :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)

        real(rk),   allocatable    :: vals(:)

        real(rk) :: gam

        !& DEBUG, HARDCODED GAMMA
        gam = 1.4_rk

        vals = (gam-ONE)*(rhoE - HALF*rho*((rhou*rhou)/(rho*rho) + (rhov*rhov)/(rho*rho) + (rhow*rhow)/(rho*rho)))

    end function compute_pressure_real
    !************************************************************************************













    !> Returns a constant gamma value
    !!
    !!
    !!   @author Nathan A. Wukie
    !!   @date   2/1/2016
    !!
    !------------------------------------------------------------------------
    function compute_gamma_ad(self,rho,rhou,rhov,rhow,rhoE) result(vals)
        class(perfect_gas_t),       intent(in)      :: self
        type(AD_D),                 intent(in)      :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)

        type(AD_D), allocatable  :: vals(:)

        !
        ! Make sure vals derivatives are initialized
        !
        vals = ZERO*rho

        !
        ! Set constant value
        !
        vals = 1.4_rk

    end function compute_gamma_ad
    !***********************************************************************




    !> Returns a constant gamma value
    !!
    !!
    !!   @author Nathan A. Wukie
    !!   @date   2/1/2016
    !!
    !------------------------------------------------------------------------
    function compute_gamma_real(self,rho,rhou,rhov,rhow,rhoE) result(vals)
        class(perfect_gas_t),       intent(in)      :: self
        real(rk),                   intent(in)      :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)

        real(rk),   allocatable :: vals(:)

        !
        ! Make sure array is allocated
        !
        vals = ZERO*rho

        !
        ! Set constant value
        !
        vals = 1.4_rk

    end function compute_gamma_real
    !***********************************************************************






    !> Returns a constant gamma value
    !!
    !!
    !!   @author Nathan A. Wukie (AFRL)
    !!   @date   9/23/2016
    !!
    !------------------------------------------------------------------------
    function compute_temperature_ad(self,rho,rhou,rhov,rhow,rhoE) result(vals)
        class(perfect_gas_t),       intent(in)      :: self
        type(AD_D),                 intent(in)      :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)

        type(AD_D), allocatable  :: vals(:), P(:)

        !
        ! Make sure vals derivatives are initialized
        !
        P = self%compute_pressure(rho,rhou,rhov,rhow,rhoE)

        !
        ! Set constant value
        !
        vals = P/(rho*self%R)

    end function compute_temperature_ad
    !***********************************************************************






    !>  Sutherlands Law for Viscosity
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/23/2016
    !!
    !------------------------------------------------------------------------
    function compute_viscosity_dynamic_ad(self,T) result(vals)
        class(perfect_gas_t),       intent(in)      :: self
        type(AD_D),                 intent(in)      :: T(:)

        type(AD_D), allocatable  :: vals(:)

        real(rk) :: mu0 = 1.7894e-5_rk  ! [kg/(m*s)]
        real(rk) :: T0  = 273.11_rk     ! [K]
        real(rk) :: S   = 110.56_rk     ! [K]

        !
        ! Set constant value
        !
        !vals = mu0*((T/T0)**(THREE/TWO))*(T0+S)/(T+S)

        vals = T
        vals = 1._rk

    end function compute_viscosity_dynamic_ad
    !***********************************************************************







    !>  Stokes' Hypothesis for Second Coefficient of Viscosity
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/23/2016
    !!
    !------------------------------------------------------------------------
    function compute_viscosity_second_ad(self,mu,T) result(vals)
        class(perfect_gas_t),       intent(in)      :: self
        type(AD_D),                 intent(in)      :: mu(:)
        type(AD_D),                 intent(in)      :: T(:)

        type(AD_D), allocatable  :: vals(:)

        !
        ! Set constant value
        !
        vals = -(TWO/THREE)*mu

    end function compute_viscosity_second_ad
    !***********************************************************************





end module perfect_gas
