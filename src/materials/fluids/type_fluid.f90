module type_fluid
    use mod_kinds,      only: rk
    use type_material,  only: material_t
    use DNAD_D
    implicit none
    private





    !>  Abstract fluid class
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/15/2016
    !!
    !--------------------------------------------------------------------------
    type, public, extends(material_t), abstract :: fluid_t


    contains

        generic, public                     :: compute_pressure          => compute_pressure_ad, compute_pressure_real
        generic, public                     :: compute_gamma             => compute_gamma_ad,    compute_gamma_real
        generic, public                     :: compute_temperature       => compute_temperature_ad
        generic, public                     :: compute_viscosity_dynamic => compute_viscosity_dynamic_ad
        generic, public                     :: compute_viscosity_second  => compute_viscosity_second_ad

        procedure(compute_ad),   deferred   :: compute_pressure_ad
        procedure(compute_real), deferred   :: compute_pressure_real

        procedure(compute_ad),   deferred   :: compute_gamma_ad
        procedure(compute_real), deferred   :: compute_gamma_real

        procedure(compute_ad),   deferred   :: compute_temperature_ad

        procedure(compute_viscosity_ad),          deferred   :: compute_viscosity_dynamic_ad
        procedure(compute_second_viscosity_ad),   deferred   :: compute_viscosity_second_ad

    end type fluid_t
    !**************************************************************************




    abstract interface
        function compute_ad(self,rho,rhou,rhov,rhow,rhoE) result(vals)
            import fluid_t
            import AD_D

            class(fluid_t),             intent(in)      :: self
            type(AD_D),                 intent(in)      :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)
            type(AD_D), allocatable :: vals(:)
        end function
    end interface



    abstract interface
        function compute_viscosity_ad(self,T) result(vals)
            import fluid_t
            import AD_D

            class(fluid_t),             intent(in)      :: self
            type(AD_D),                 intent(in)      :: T(:)
            type(AD_D), allocatable :: vals(:)
        end function
    end interface


    abstract interface
        function compute_second_viscosity_ad(self,mu,T) result(vals)
            import fluid_t
            import AD_D

            class(fluid_t),             intent(in)      :: self
            type(AD_D),                 intent(in)      :: mu(:)
            type(AD_D),                 intent(in)      :: T(:)
            type(AD_D), allocatable :: vals(:)
        end function
    end interface




    abstract interface
        function compute_real(self,rho,rhou,rhov,rhow,rhoE) result(vals)
            use mod_kinds,  only: rk
            import fluid_t

            class(fluid_t),             intent(in)      :: self
            real(rk),                   intent(in)      :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)
            real(rk),   allocatable :: vals(:)
        end function
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






























end module type_fluid
