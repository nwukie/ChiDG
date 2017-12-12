module model_pressure_gradient
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE
    use mod_fluid,          only: Rgas, gam
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !>  An equation of state model for an ideal gas.
    !!
    !!  Model Fields:
    !!      - Pressure Gradient - 1
    !!      - Pressure Gradient - 2
    !!      - Pressure Gradient - 3
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/6/2017
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: pressure_gradient_t

    contains

        procedure   :: init
        procedure   :: compute

    end type pressure_gradient_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/6/2017
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(pressure_gradient_t), intent(inout)   :: self

        call self%set_name('Pressure Gradient')
        !call self%set_dependency('f(Q-)')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Pressure Gradient - 1')
        call self%add_model_field('Pressure Gradient - 2')
        call self%add_model_field('Pressure Gradient - 3')


    end subroutine init
    !***************************************************************************************






    !>  Routine for computing pressure gradient.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/6/2017
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(pressure_gradient_t), intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable ::                              &
            density, mom1, mom2, mom3, energy,                                  &
            grad1_density, grad1_mom1, grad1_mom2, grad1_mom3, grad1_energy,    &
            grad2_density, grad2_mom1, grad2_mom2, grad2_mom3, grad2_energy,    &
            grad3_density, grad3_mom1, grad3_mom2, grad3_mom3, grad3_energy,    &
            dp_ddensity, dp_dmom1, dp_dmom2, dp_dmom3, dp_denergy,              &
            grad1_p, grad2_p, grad3_p

        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_field('Density',    'value')
        mom1    = worker%get_field('Momentum-1', 'value')
        mom2    = worker%get_field('Momentum-2', 'value')
        mom3    = worker%get_field('Momentum-3', 'value')
        energy  = worker%get_field('Energy',     'value')

        grad1_density = worker%get_field('Density',    'grad1')
        grad1_mom1    = worker%get_field('Momentum-1', 'grad1')
        grad1_mom2    = worker%get_field('Momentum-2', 'grad1')
        grad1_mom3    = worker%get_field('Momentum-3', 'grad1')
        grad1_energy  = worker%get_field('Energy',     'grad1')


        grad2_density = worker%get_field('Density',    'grad2')
        grad2_mom1    = worker%get_field('Momentum-1', 'grad2')
        grad2_mom2    = worker%get_field('Momentum-2', 'grad2')
        grad2_mom3    = worker%get_field('Momentum-3', 'grad2')
        grad2_energy  = worker%get_field('Energy',     'grad2')


        grad3_density = worker%get_field('Density',    'grad3')
        grad3_mom1    = worker%get_field('Momentum-1', 'grad3')
        grad3_mom2    = worker%get_field('Momentum-2', 'grad3')
        grad3_mom3    = worker%get_field('Momentum-3', 'grad3')
        grad3_energy  = worker%get_field('Energy',     'grad3')



        dp_ddensity =  (gam-ONE)*HALF*(mom1*mom1 + mom2*mom2 + mom3*mom3)/(density*density)
        dp_dmom1    = -(gam-ONE)*mom1/density
        dp_dmom2    = -(gam-ONE)*mom2/density
        dp_dmom3    = -(gam-ONE)*mom3/density
        dp_denergy  = dp_ddensity ! init storage
        dp_denergy  =  (gam-ONE)


        ! Compute pressure gradient using Chain-rule
        grad1_p = dp_ddensity * grad1_density  + &
                  dp_dmom1    * grad1_mom1     + &
                  dp_dmom2    * grad1_mom2     + &
                  dp_dmom3    * grad1_mom3     + &
                  dp_denergy  * grad1_energy

        grad2_p = dp_ddensity * grad2_density  + &
                  dp_dmom1    * grad2_mom1     + &
                  dp_dmom2    * grad2_mom2     + &
                  dp_dmom3    * grad2_mom3     + &
                  dp_denergy  * grad2_energy

        grad3_p = dp_ddensity * grad3_density  + &
                  dp_dmom1    * grad3_mom1     + &
                  dp_dmom2    * grad3_mom2     + &
                  dp_dmom3    * grad3_mom3     + &
                  dp_denergy  * grad3_energy



        call worker%store_model_field('Pressure Gradient - 1', 'value', grad1_p)
        call worker%store_model_field('Pressure Gradient - 2', 'value', grad2_p)
        call worker%store_model_field('Pressure Gradient - 3', 'value', grad3_p)


    end subroutine compute
    !***************************************************************************************









end module model_pressure_gradient
