module euler_ale_volume_operator
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF

    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    use ieee_arithmetic,        only: ieee_is_nan
    implicit none

    private

    
    !> Volume flux for Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: euler_ale_volume_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type euler_ale_volume_operator_t
    !******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_ale_volume_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("Euler ALE Volume Flux")

        ! Set operator type
        call self%set_operator_type("Volume Advective Flux")

        ! Set operator equations
        call self%add_primary_field("Density"   )
        call self%add_primary_field("Momentum-1")
        call self%add_primary_field("Momentum-2")
        call self%add_primary_field("Momentum-3")
        call self%add_primary_field("Energy"    )

    end subroutine init
    !********************************************************************************



    !> Volume flux routine for Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(euler_ale_volume_operator_t), intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:) ::            &
            density, mom1, mom2, mom3, energy, p, enthalpy, &
            flux_1, flux_2, flux_3, invdensity


        type(AD_D), allocatable, dimension(:,:)  :: flux


        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_primary_field_element('Density'   , 'value')
        mom1    = worker%get_primary_field_element('Momentum-1', 'value')
        mom2    = worker%get_primary_field_element('Momentum-2', 'value')
        mom3    = worker%get_primary_field_element('Momentum-3', 'value')
        energy  = worker%get_primary_field_element('Energy'    , 'value')


        invdensity = ONE/density
    


        !
        ! Compute pressure and total enthalpy
        !
        p = worker%get_model_field_element('Pressure','value')

        enthalpy = (energy + p)*invdensity

        !===========================
        !        MASS FLUX
        !===========================
        flux_1 = mom1
        flux_2 = mom2
        flux_3 = mom3

        call worker%integrate_volume_flux_ale('Density','Advective',flux_1,flux_2,flux_3)

        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux_1 = (mom1*mom1)*invdensity  +  p
        flux_2 = (mom1*mom2)*invdensity
        flux_3 = (mom1*mom3)*invdensity
        
        call worker%integrate_volume_flux_ale('Momentum-1','Advective',flux_1,flux_2,flux_3)

        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux_1 = (mom2*mom1)*invdensity
        flux_2 = (mom2*mom2)*invdensity  +  p
        flux_3 = (mom2*mom3)*invdensity
        
        call worker%integrate_volume_flux_ale('Momentum-2','Advective',flux_1,flux_2,flux_3)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        flux_1 = (mom3*mom1)*invdensity
        flux_2 = (mom3*mom2)*invdensity
        flux_3 = (mom3*mom3)*invdensity  +  p

        call worker%integrate_volume_flux_ale('Momentum-3','Advective',flux_1,flux_2,flux_3)

        !============================
        !       ENERGY FLUX
        !============================
        flux_1 = mom1*enthalpy
        flux_2 = mom2*enthalpy
        flux_3 = mom3*enthalpy

        call worker%integrate_volume_flux_ale('Energy','Advective',flux_1,flux_2,flux_3)

    end subroutine compute
    !*********************************************************************************************************






end module euler_ale_volume_operator
