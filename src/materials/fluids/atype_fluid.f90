module atype_fluid
    use DNAD_D
    implicit none
    private

    !
    !
    !
    !
    !---------------------------------------------------------------
    type, public, abstract :: fluid_t


    contains
        procedure(compute), deferred :: compute_pressure
        procedure(compute), deferred :: compute_gamma
    end type fluid_t
    !---------------------------------------------------------------




    abstract interface
        subroutine compute(self,rho,rhou,rhov,rhow,rhoE,vals)
            import fluid_t
            import AD_D

            class(fluid_t), intent(in)      :: self
            type(AD_D),     intent(in)      :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)
            type(AD_D),     intent(inout)   :: vals(:)
        end subroutine
    end interface




end module atype_fluid
