module perfect_gas
    use mod_kinds,      only: rk
    use mod_constants,  only: ONE, HALF
    use atype_fluid,    only: fluid_t
    use DNAD_D


    ! A fluid type under perfect gas assumptions.
    !
    !   @author Nathan A. Wukie
    !-----------------------------------------------------------------------
    type, extends(fluid_t), public :: perfect_gas_t

    contains
        procedure :: compute_pressure

    end type perfect_gas_t


contains



    !   Compute pressure of the fluid, under the perfect gas assumption
    !
    !   @author Nathan A. Wukie
    !
    !   @param[in]      rho     Fluid density
    !   @param[in]      rhou    Fluid x-momentum
    !   @param[in]      rhov    Fluid y-momentum
    !   @param[in]      rhow    Fluid z-momentum
    !   @param[in]      rhoE    Fluid Total Energy
    !   @param[inout]   p       Fluid pressure
    !------------------------------------------------------------------------
    subroutine compute_pressure(self,rho,rhou,rhov,rhow,rhoE,p)
        class(perfect_gas_t),   intent(in)  :: self
        type(AD_D),             intent(in)      :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)
        type(AD_D),             intent(inout)   :: p(:)

        real(rk)    :: gam

        gam = 1.4_rk

        p = (gam-ONE)*(rhoE - HALF*rho*((rhou*rhou)/(rho*rho) + (rhov*rhov)/(rho*rho) + (rhow*rhow)/(rho*rho)))

    end subroutine
    !------------------------------------------------------------------------

end module
