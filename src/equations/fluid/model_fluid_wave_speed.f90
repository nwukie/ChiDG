module model_fluid_wave_speed
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, TWO
    use mod_fluid,          only: Rgas, gam
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
        call self%set_dependency('f(Q-)')

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
            density, energy, T, P, vmag, sound_speed, wave_speed


        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_field('Density', 'value')
        energy  = worker%get_field('Energy',  'value')


        !
        ! Get Temperature defined by some model for equation of state
        !
        T = worker%get_field('Temperature', 'value')
        P = worker%get_field('Pressure', 'value')



        !
        ! Compute acoustic and convective wave speeds
        !   Danger computing vmag with sqrt because derivative
        !   is undefined at zero.
        !
        !vmag = sqrt((densityu**TWO + densityv**TWO + densityw**TWO)/(density*density))
        vmag = energy - (P/(gam - ONE))
        sound_speed = sqrt(gam * Rgas * T)

        
        !
        ! Compute maximum wave speed
        !
        wave_speed = vmag + sound_speed

        call worker%store_model_field('Maximum Wave Speed', 'value', wave_speed)


    end subroutine compute
    !***************************************************************************************




end module model_fluid_wave_speed
