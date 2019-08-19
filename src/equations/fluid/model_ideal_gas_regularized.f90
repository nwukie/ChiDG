module model_ideal_gas_regularized
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, ZERO
    use mod_fluid,          only: Rgas, gam
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !> Regularized ideal gas model that enforces positivity of pressure and temperature.
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    08/05/2019 
    !!
    !--------------------------------------------------------------------------------
    type, extends(model_t)  :: ideal_gas_regularized_t

        real(rk)    :: window_size_pressure      = 1.0e-2_rk
        real(rk)    :: minimum_value_pressure    = 1.0e-11_rk
        logical     :: regularize_pressure       = .true.


    contains

        procedure   :: init
        procedure   :: compute

    end type ideal_gas_regularized_t
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
        class(ideal_gas_regularized_t), intent(inout)   :: self

        integer             :: unit, msg
        logical             :: file_exists
        real(rk)    :: window_size_pressure, minimum_value_pressure
        logical     :: regularize_pressure

        namelist /regularize_ideal_gas/ window_size_pressure, minimum_value_pressure, regularize_pressure

        call self%set_name('Regularized Ideal Gas')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('Pressure')
        call self%add_model_field('Temperature')

        inquire(file='regularize_ideal_gas.nml', exist=file_exists)
        if (file_exists) then
            open(newunit=unit,form='formatted',file='regularize_ideal_gas.nml')
            read(unit,nml=regularize_ideal_gas,iostat=msg)
            if (msg == 0) self%window_size_pressure = window_size_pressure
            if (msg == 0) self%minimum_value_pressure = minimum_value_pressure
            if (msg == 0) self%regularize_pressure = regularize_pressure
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
        class(ideal_gas_regularized_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            density, mom1, mom2, mom3, energy, pressure, temperature


        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_field('Density',    'value')
        mom1    = worker%get_field('Momentum-1', 'value')
        mom2    = worker%get_field('Momentum-2', 'value')
        mom3    = worker%get_field('Momentum-3', 'value')
        energy  = worker%get_field('Energy',     'value')


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2 = mom2 / worker%coordinate('1')
        end if


        pressure = (gam-ONE)*(energy - HALF*( (mom1*mom1) + (mom2*mom2) + (mom3*mom3) )/density)
        if (self%regularize_pressure) pressure = pressure*sin_ramp(pressure, ZERO, self%window_size_pressure) + self%minimum_value_pressure
        temperature = pressure/(density*Rgas)



        call worker%store_model_field('Pressure',    'value', pressure)
        call worker%store_model_field('Temperature', 'value', temperature)



    end subroutine compute
    !***************************************************************************************









end module model_ideal_gas_regularized
