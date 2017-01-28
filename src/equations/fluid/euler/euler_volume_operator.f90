module euler_volume_operator
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF

    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    implicit none

    private

    
    !> Volume flux for Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: euler_volume_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type euler_volume_operator_t
    !******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_volume_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("Euler Volume Flux")

        ! Set operator type
        call self%set_operator_type("Volume Advective Flux")

        ! Set operator equations
        call self%add_primary_field("Density"   )
        call self%add_primary_field("X-Momentum")
        call self%add_primary_field("Y-Momentum")
        call self%add_primary_field("Z-Momentum")
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
        class(euler_volume_operator_t), intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop



        type(AD_D), allocatable, dimension(:) ::    &
            rho, rhou, rhov, rhow, rhoE, p, H,      &
            flux_x, flux_y, flux_z, invrho



        !
        ! Interpolate solution to quadrature nodes
        !
        rho  = worker%get_primary_field_element("Density"   ,'value')
        rhou = worker%get_primary_field_element("X-Momentum",'value')
        rhov = worker%get_primary_field_element("Y-Momentum",'value')
        rhow = worker%get_primary_field_element("Z-Momentum",'value')
        rhoE = worker%get_primary_field_element("Energy"    ,'value')


        invrho = ONE/rho



        !
        ! Compute pressure and total enthalpy
        !
        p = worker%get_model_field_element('Pressure','value')
        H = (rhoE + p)*invrho

        !===========================
        !        MASS FLUX
        !===========================
        flux_x = rhou
        flux_y = rhov
        flux_z = rhow

        call worker%integrate_volume('Density',flux_x,flux_y,flux_z)


        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux_x = (rhou*rhou)*invrho  +  p
        flux_y = (rhou*rhov)*invrho
        flux_z = (rhou*rhow)*invrho

        call worker%integrate_volume('X-Momentum',flux_x,flux_y,flux_z)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux_x = (rhov*rhou)*invrho
        flux_y = (rhov*rhov)*invrho  +  p
        flux_z = (rhov*rhow)*invrho

        call worker%integrate_volume('Y-Momentum',flux_x,flux_y,flux_z)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        flux_x = (rhow*rhou)*invrho
        flux_y = (rhow*rhov)*invrho
        flux_z = (rhow*rhow)*invrho  +  p

        call worker%integrate_volume('Z-Momentum',flux_x,flux_y,flux_z)

        !============================
        !       ENERGY FLUX
        !============================
        flux_x = rhou*H
        flux_y = rhov*H
        flux_z = rhow*H

        call worker%integrate_volume('Energy',flux_x,flux_y,flux_z)

    end subroutine compute
    !*********************************************************************************************************






end module euler_volume_operator
