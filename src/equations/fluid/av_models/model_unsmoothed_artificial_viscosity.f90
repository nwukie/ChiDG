module model_unsmoothed_artificial_viscosity
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
    type, extends(model_t)  :: unsmoothed_artificial_viscosity_t


    contains

        procedure   :: init
        procedure   :: compute

    end type unsmoothed_artificial_viscosity_t
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
        class(unsmoothed_artificial_viscosity_t), intent(inout)   :: self
        




        call self%set_name('Unsmoothed Artificial Viscosity')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Smoothed Artificial Bulk Viscosity')
        call self%add_model_field('Smoothed Artificial Shear Viscosity')
        call self%add_model_field('Smoothed Artificial Thermal Conductivity')

        

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
        class(unsmoothed_artificial_viscosity_t),   intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            density, vel1, vel2, vel3, T, c, wave_speed, shock_sensor, &
            art_bulk_viscosity, art_shear_viscosity, art_thermal_conductivity

        art_bulk_viscosity          =  worker%get_field('Artificial Bulk Viscosity', 'value')
        art_shear_viscosity         =  worker%get_field('Artificial Shear Viscosity', 'value')
        art_thermal_conductivity    =  worker%get_field('Artificial Thermal Conductivity', 'value')




        !
        ! Contribute laminar viscosity
        !
        call worker%store_model_field('Smoothed Artificial Bulk Viscosity', 'value', art_bulk_viscosity)
        call worker%store_model_field('Smoothed Artificial Shear Viscosity', 'value', art_shear_viscosity)
        call worker%store_model_field('Smoothed Artificial Thermal Conductivity', 'value', art_thermal_conductivity)


    end subroutine compute
    !***************************************************************************************




end module model_unsmoothed_artificial_viscosity
