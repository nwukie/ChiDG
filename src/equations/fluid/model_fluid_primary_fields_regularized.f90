module model_fluid_primary_fields_regularized
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, ZERO
    use mod_fluid,          only: Rgas, gam
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !> This model computes regularized values for fluid primary field variables by
    !! enforcing positivity of density and energy.
    !!
    !! @author  Eric M. Wolf
    !! @date    08/05/2019 
    !!
    !--------------------------------------------------------------------------------
    type, extends(model_t)  :: fluid_primary_fields_regularized_t

        real(rk)    :: window_size_density      = 1.0e-2_rk
        real(rk)    :: minimum_value_density    = 1.0e-11_rk
        logical     :: regularize_density       = .true.

        real(rk)    :: window_size_energy       = 1.0e-2_rk
        real(rk)    :: minimum_value_energy     = 1.0e-11_rk
        logical     :: regularize_energy        = .true.

    contains

        procedure   :: init
        procedure   :: compute

    end type fluid_primary_fields_regularized_t
    !***************************************************************************************





contains




    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    08/05/2019 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)   
        class(fluid_primary_fields_regularized_t), intent(inout)   :: self

        integer             :: unit, msg
        logical             :: file_exists
        real(rk)    :: window_size_density, minimum_value_density, window_size_energy, minimum_value_energy
        logical     :: regularize_density, regularize_energy

        namelist /regularize_fluid_primary_fields/ window_size_density, minimum_value_density, regularize_density, &
            window_size_energy, minimum_value_energy, regularize_energy

        call self%set_name('Regularized Fluid Primary Fields')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('Density')
        call self%add_model_field('Energy')

        inquire(file='regularize_fluid_primary_fields.nml', exist=file_exists)
        if (file_exists) then
            open(newunit=unit,form='formatted',file='regularize_fluid_primary_fields.nml')
            read(unit,nml=regularize_fluid_primary_fields,iostat=msg)
            if (msg == 0) self%window_size_density = window_size_density
            if (msg == 0) self%minimum_value_density = minimum_value_density
            if (msg == 0) self%regularize_density = regularize_density
            if (msg == 0) self%window_size_energy = window_size_energy
            if (msg == 0) self%minimum_value_energy = minimum_value_energy
            if (msg == 0) self%regularize_energy = regularize_energy
            close(unit)
        end if

    end subroutine init
    !***************************************************************************************




    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    08/05/2019 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(fluid_primary_fields_regularized_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            density, energy


        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_field('Density',    'value')
        energy  = worker%get_field('Energy',     'value')

        if (self%regularize_density) density = density*sin_ramp(density, ZERO, self%window_size_density) + self%minimum_value_density
        if (self%regularize_energy) energy = energy*sin_ramp(energy, ZERO, self%window_size_energy) + self%minimum_value_energy



        call worker%store_model_field('Density',    'value', density)
        call worker%store_model_field('Energy', 'value', energy)



    end subroutine compute
    !***************************************************************************************









end module model_fluid_primary_fields_regularized
