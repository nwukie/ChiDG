module model_zero_artificial_viscosity
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
    type, extends(model_t)  :: zero_artificial_viscosity_t

        real(rk) :: av_constant = 1.5_rk

    contains

        procedure   :: init
        procedure   :: compute

    end type zero_artificial_viscosity_t
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
        class(zero_artificial_viscosity_t), intent(inout)   :: self
        
        real(rk)            :: av_constant
        integer             :: unit, msg
        logical             :: file_exists

        namelist /zero_artificial_viscosity/   av_constant



        call self%set_name('Zero Artificial Viscosity')
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
            read(unit,nml=zero_artificial_viscosity,iostat=msg)
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
        class(zero_artificial_viscosity_t),   intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            density, vel1, vel2, vel3, T, c, wave_speed, shock_sensor, &
            art_bulk_viscosity, art_shear_viscosity, art_thermal_conductivity

        real(rk), dimension(:)  :: h(3)
        real(rk)                :: hmin
        real(rk)                :: Pr_star  = 0.9_rk     
        integer(ik)             :: p, ii
        

                
        density = worker%get_field('Density','value')

        art_bulk_viscosity = ZERO*density

        art_shear_viscosity = ZERO*density

        art_thermal_conductivity = ZERO*density

        

        !
        ! Contribute laminar viscosity
        !
        call worker%store_model_field('Artificial Bulk Viscosity', 'value', art_bulk_viscosity)
        call worker%store_model_field('Artificial Shear Viscosity', 'value', art_shear_viscosity)
        call worker%store_model_field('Artificial Thermal Conductivity', 'value', art_thermal_conductivity)


    end subroutine compute
    !***************************************************************************************




end module model_zero_artificial_viscosity
