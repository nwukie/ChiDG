module type_fluid_wave_speed
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, TWO
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !>  A model defining the maximum wave speed of the fluid system of equations.
    !!
    !!  Model Fields:
    !!      - Maximum Wave Speed
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: fluid_wave_speed_t

        real(rk)    :: gam = 1.4_rk     ! ratio of specific heats
        real(rk)    :: R   = 287.15_rk  ! ideal gas constant [J/(kg*K)]

    contains

        procedure   :: init
        procedure   :: compute

    end type fluid_wave_speed_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/31/2017
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(fluid_wave_speed_t), intent(inout)   :: self

        call self%set_name('Fluid Wave Speed')
        call self%set_dependency('Q-')

        call self%add_model_field('Maximum Wave Speed')

    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Pressure and Temperature.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(fluid_wave_speed_t),    intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            rho, rhoE, T, P, vmag, sound_speed, wave_speed


        !
        ! Interpolate solution to quadrature nodes
        !
        rho  = worker%get_primary_field_general('Density',    'value')
        rhoE = worker%get_primary_field_general('Energy', 'value')


        !
        ! Get Temperature defined by some model for equation of state
        !
        T    = worker%get_model_field_general('Temperature', 'value')
        P    = worker%get_model_field_general('Pressure', 'value')



        !
        ! Compute acoustic and convective wave speeds
        !   Danger computing vmag with sqrt because derivative
        !   is undefined at zero.
        !
        !vmag = sqrt((rhou**TWO + rhov**TWO + rhow**TWO)/(rho*rho))
        vmag = rhoE - (P/(self%gam - ONE))
        sound_speed = sqrt(self%gam * self%R * T)

        
        !
        ! Compute maximum wave speed
        !
        wave_speed = vmag + sound_speed

        call worker%store_model_field('Maximum Wave Speed', 'value', wave_speed)


    end subroutine compute
    !***************************************************************************************




end module type_fluid_wave_speed
