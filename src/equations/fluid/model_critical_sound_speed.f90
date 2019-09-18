module model_critical_sound_speed
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: ZERO, HALF, ONE, TWO, PI, RKTOL
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use mod_fluid
    use DNAD_D
    use ieee_arithmetic 
    implicit none


    


    !> Based on the reference:
    !! Fernandez, Pablo, Cuong Nguyen, and Jaime Peraire. 
    !! "A physics-based shock capturing method for unsteady laminar and turbulent flows." 
    !! 2018 AIAA Aerospace Sciences Meeting. 2018.
    !!
    !! @author  Eric M. Wolf
    !! @date    09/06/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(model_t)  :: critical_sound_speed_t


    contains

        procedure   :: init
        procedure   :: compute

    end type critical_sound_speed_t
    !***************************************************************************************





contains




    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    09/06/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)   
        class(critical_sound_speed_t), intent(inout)   :: self

        call self%set_name('Critical Sound Speed')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('Critical Sound Speed')
        call self%add_model_field('Stagnation Temperature')
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
        class(critical_sound_speed_t),  intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        integer(ik) :: ii
        type(AD_D), allocatable,    dimension(:) ::         &
            t, t_stag, c_star,   &
            div_vel, density, mom1, mom2, mom3, vel_mag

        real(rk),   allocatable,    dimension(:) :: h
        real(rk) :: Dmax, d0, delta_d

        density = worker%get_field('Density',    'value')
        mom1    = worker%get_field('Momentum-1', 'value')
        mom2    = worker%get_field('Momentum-2', 'value')
        mom3    = worker%get_field('Momentum-3', 'value')

        vel_mag = ZERO*density
        vel_mag = (mom1**TWO+mom2**TWO+mom3**TWO)/density**TWO 
        t = worker%get_field('Temperature', 'value') 
        t_stag = t + vel_mag/(TWO*cp)

        call worker%store_model_field('Stagnation Temperature', 'value', t_stag)

        c_star = ZERO*density
        c_star = sqrt(gam*Rgas*(TWO/(gam+ONE))*((t_stag)))

        call worker%store_model_field('Critical Sound Speed', 'value', c_star)

    end subroutine compute
    !***************************************************************************************




end module model_critical_sound_speed
