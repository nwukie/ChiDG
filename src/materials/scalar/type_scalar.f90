module type_scalar
    use mod_kinds,      only: rk
    use type_material,  only: material_t
    use DNAD_D
    implicit none
    private





    !>  Abstract scalar class
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/15/2016
    !!
    !--------------------------------------------------------------------------
    type, public, extends(material_t), abstract :: scalar_t


    contains

        procedure(compute_ad),   deferred   :: compute_mu

    end type scalar_t
    !**************************************************************************




    abstract interface
        subroutine compute_ad(self,u,vals)
            import scalar_t
            import AD_D

            class(scalar_t),            intent(in)      :: self
            type(AD_D),                 intent(in)      :: u(:)
            type(AD_D), allocatable,    intent(inout)   :: vals(:)
        end subroutine
    end interface



    abstract interface
        subroutine compute_real(self,rho,rhou,rhov,rhow,rhoE,vals)
            use mod_kinds,  only: rk
            import scalar_t

            class(scalar_t),             intent(in)      :: self
            real(rk),                   intent(in)      :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)
            real(rk),   allocatable,    intent(inout)   :: vals(:)
        end subroutine
    end interface

contains




!    !>
!    !!
!    !!  @author Nathan A. Wukie (AFRL)
!    !!  @date   8/30/2016
!    !!
!    !!
!    !------------------------------------------------------------------------------------------
!    subroutine init(self)
!        class(fluid_t), intent(inout)   :: self
!
!
!        call self%contribute_equation("Density")
!        call self%contribute_equation("X-Momentum")
!        call self%contribute_equation("Y-Momentum")
!        call self%contribute_equation("Z-Momentum")
!        call self%contribute_equation("Energy")
!        call self%contribute_equation("Turbulence Viscosity")
!        call self%contribute_equation("Turbulence Kinetic Energy")
!        call self%contribute_equation("Turbulence Dissipation Rate")
!
!
!    end subroutine init
!    !******************************************************************************************






























end module type_scalar
