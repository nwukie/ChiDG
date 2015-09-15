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
        procedure(pressure), deferred :: compute_pressure
    end type fluid_t
    !---------------------------------------------------------------




    abstract interface
        subroutine pressure(self,rho,rhou,rhov,rhow,rhoE,p)
            import fluid_t
            import AD_D

            class(fluid_t), intent(in)      :: self
            type(AD_D),     intent(in)      :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)
            type(AD_D),     intent(inout)   :: p(:)
        end subroutine
    end interface




end module atype_fluid
