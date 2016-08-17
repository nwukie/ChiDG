module type_fluid
    use mod_kinds,  only: rk
    use DNAD_D
    implicit none
    private





    !>  Abstract fluid class
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/15/2016
    !!
    !--------------------------------------------------------------------------
    type, public, abstract :: fluid_t


    contains

        generic, public                     :: compute_pressure => compute_pressure_ad, compute_pressure_real
        generic, public                     :: compute_gamma    => compute_gamma_ad,    compute_gamma_real

        procedure(compute_ad),   deferred   :: compute_pressure_ad
        procedure(compute_real), deferred   :: compute_pressure_real

        procedure(compute_ad),   deferred   :: compute_gamma_ad
        procedure(compute_real), deferred   :: compute_gamma_real

    end type fluid_t
    !**************************************************************************




    abstract interface
        subroutine compute_ad(self,rho,rhou,rhov,rhow,rhoE,vals)
            import fluid_t
            import AD_D

            class(fluid_t), intent(in)      :: self
            type(AD_D),     intent(in)      :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)
            type(AD_D),     intent(inout)   :: vals(:)
        end subroutine
    end interface



    abstract interface
        subroutine compute_real(self,rho,rhou,rhov,rhow,rhoE,vals)
            use mod_kinds,  only: rk
            import fluid_t

            class(fluid_t), intent(in)      :: self
            real(rk),       intent(in)      :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)
            real(rk),       intent(inout)   :: vals(:)
        end subroutine
    end interface


end module type_fluid
