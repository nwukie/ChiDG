module model_artificial_viscosity
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: THREE, TWO, ONE, ZERO
    use mod_fluid
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    use ieee_arithmetic
    implicit none


    


    !>  Artificial viscosity model based on
    !!  A physics-based shock capturing method for unsteady laminar and turbulent flows
    !!      Fernandez et al, 2018, AIAA SciTech Forum
    !!
    !!  Model Fields:
    !!      - Artifical Bulk Viscosity
    !!      - Artifical Shear Viscosity
    !!      - Artifical Thermal Conductivity 
    !!
    !!  @author Eric M. Wolf
    !!  @date   07/11/2018
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: artificial_viscosity_t

        real(rk) :: av_constant = 1.5_rk

    contains

        procedure   :: init
        procedure   :: compute

    end type artificial_viscosity_t
    !***************************************************************************************





contains




    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    07/11/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)   
        class(artificial_viscosity_t), intent(inout)   :: self
        
        real(rk)            :: av_constant
        integer             :: unit, msg
        logical             :: file_exists

        namelist /artificial_viscosity/   av_constant



        call self%set_name('Artificial Viscosity')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Artificial Bulk Viscosity')
        call self%add_model_field('Artificial Shear Viscosity')
        call self%add_model_field('Artificial Thermal Conductivity')

        !
        ! Check if input from 'models.nml' is available.
        !   1: if available, read and set self%mu
        !   2: if not available, do nothing and mu retains default value
        !
        inquire(file='models.nml', exist=file_exists)
        if (file_exists) then
            open(newunit=unit,form='formatted',file='models.nml')
            read(unit,nml=artificial_viscosity,iostat=msg)
            if (msg == 0) self%av_constant = av_constant
            close(unit)
        end if


    end subroutine init
    !***************************************************************************************





    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    07/11/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(artificial_viscosity_t),   intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            density, vel1, vel2, vel3, T, c, wave_speed, shock_sensor, &
            art_bulk_viscosity, art_shear_viscosity, art_thermal_conductivity

        real(rk), dimension(:)  :: h(3)
        real(rk)                :: hmin
        real(rk)                :: Pr_star  = 0.9_rk     
        integer(ik)             :: p, ii
        

        h = worker%element_size('interior')
        p = worker%solution_order('interior')
        hmin = minval(h)/real((p+ONE), rk)
        
        
        density = worker%get_field('Density','value')

        vel1 = worker%get_field('Momentum-1','value')/density
        vel2 = worker%get_field('Momentum-2','value')/density
        vel3 = worker%get_field('Momentum-3','value')/density

        T = worker%get_field('Temperature','value')
        if (any(ieee_is_nan(T(:)%x_ad_))) print *, 'av: temp is nan, interp source: ', worker%interpolation_source
        c = T
        c = sqrt((gam*Rgas*T))
        if (any(ieee_is_nan(c(:)%x_ad_))) print *, 'av: c is nan, interp source: ', worker%interpolation_source
        if (any(ieee_is_nan(c(:)%x_ad_))) print *, 'temp: ', T(:)%x_ad_

        wave_speed = c
        wave_speed = sqrt(vel1**TWO+vel2**TWO+vel3**TWO+c**TWO)

        shock_sensor = worker%get_field('Shock Sensor','value')
        if (any(ieee_is_nan(density(:)%x_ad_))) print *, 'av: density is nan, interp source: ', worker%interpolation_source
        if (any(ieee_is_nan(wave_speed(:)%x_ad_))) print *, 'av: wave_speed is nan, interp source: ', worker%interpolation_source
        art_bulk_viscosity = density*self%av_constant*hmin*wave_speed*shock_sensor

        art_shear_viscosity = art_bulk_viscosity

        art_thermal_conductivity = (cp/Pr_star*art_bulk_viscosity + cp*art_bulk_viscosity)

        
        if (any(ieee_is_nan(shock_sensor(:)%x_ad_))) print *, 'shock sensor model is nan, interp source: ', worker%interpolation_source
        if (any(ieee_is_nan(art_bulk_viscosity(:)%x_ad_))) print *, 'art bulk visc model is nan'
        if (any(ieee_is_nan(art_shear_viscosity(:)%x_ad_))) print *, 'art shear visc model is nan'
        if (any(ieee_is_nan(art_thermal_conductivity(:)%x_ad_))) print *, 'art thermal conductivity model  is nan'
        do ii = 1, size(art_bulk_viscosity)
            
            !if (any(ieee_is_nan(art_bulk_viscosity(ii)%xp_ad_(:)))) print *, 'art bulk visc deriv is nan'
            !if (any(ieee_is_nan(art_shear_viscosity(ii)%xp_ad_(:)))) print *, 'art shear visc deriv is nan'
            !if (any(ieee_is_nan(art_thermal_conductivity(ii)%xp_ad_(:)))) print *, 'art thermal conductivity deriv is nan'
            art_bulk_viscosity(ii)%xp_ad_ = 0.0_rk
            art_shear_viscosity(ii)%xp_ad_ = 0.0_rk
            art_thermal_conductivity(ii)%xp_ad_ = 0.0_rk
        end do

        !
        ! Contribute laminar viscosity
        !
        call worker%store_model_field('Artificial Bulk Viscosity', 'value', art_bulk_viscosity)
        call worker%store_model_field('Artificial Shear Viscosity', 'value', art_shear_viscosity)
        call worker%store_model_field('Artificial Thermal Conductivity', 'value', art_thermal_conductivity)


    end subroutine compute
    !***************************************************************************************




end module model_artificial_viscosity
