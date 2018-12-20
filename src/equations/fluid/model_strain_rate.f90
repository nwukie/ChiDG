module model_strain_rate
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, TWO
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !>  A model for computing shear stress
    !!
    !!  Model Fields:
    !!      : shear_11, shear_12, shear_13
    !!      :           shear_22, shear_23
    !!      :                     shear_33
    !!
    !!  Lower-triangular components are not computed because the tensor is symmetric.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: strain_rate_t


    contains

        procedure   :: init
        procedure   :: compute

    end type strain_rate_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(strain_rate_t), intent(inout)   :: self

        call self%set_name('Strain Rate')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Strain Rate-11')
        call self%add_model_field('Strain Rate-22')
        call self%add_model_field('Strain Rate-33')
        call self%add_model_field('Strain Rate-12')
        call self%add_model_field('Strain Rate-13')
        call self%add_model_field('Strain Rate-23')
        call self%add_model_field('Strain Rate-Trace')

    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Pressure and Temperature.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(strain_rate_t),  intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), allocatable,    dimension(:) ::         &
            density, mom1, mom2, mom3,                      &
            grad1_density, grad2_density, grad3_density,    &
            grad1_vel1,    grad2_vel1,    grad3_vel1,       &
            grad1_vel2,    grad2_vel2,    grad3_vel2,       &
            grad1_vel3,    grad2_vel3,    grad3_vel3,       &
            strain_rate_11, strain_rate_22, strain_rate_33, &
            strain_rate_12, strain_rate_13, strain_rate_23, &
            strain_rate_trace

        real(rk),   allocatable,    dimension(:) :: r


        grad1_vel1       = worker%get_field('Velocity 1 - Gradient 1', 'value')
        grad2_vel1       = worker%get_field('Velocity 1 - Gradient 2', 'value')
        grad3_vel1       = worker%get_field('Velocity 1 - Gradient 3', 'value')

        grad1_vel2       = worker%get_field('Velocity 2 - Gradient 1', 'value')
        grad2_vel2       = worker%get_field('Velocity 2 - Gradient 2', 'value')
        grad3_vel2       = worker%get_field('Velocity 2 - Gradient 3', 'value')

        grad1_vel3       = worker%get_field('Velocity 3 - Gradient 1', 'value')
        grad2_vel3       = worker%get_field('Velocity 3 - Gradient 2', 'value')
        grad3_vel3       = worker%get_field('Velocity 3 - Gradient 3', 'value')

        strain_rate_11 = grad1_vel1
        strain_rate_22 = grad2_vel2
        strain_rate_33 = grad3_vel3

        strain_rate_12 = HALF*(grad1_vel2 + grad2_vel1)
        strain_rate_13 = HALF*(grad1_vel3 + grad3_vel1)
        strain_rate_23 = HALF*(grad2_vel3 + grad3_vel2)

        strain_rate_trace = strain_rate_11 + strain_rate_22 + strain_rate_33

        call worker%store_model_field('Strain Rate-11', 'value', strain_rate_11)
        call worker%store_model_field('Strain Rate-22', 'value', strain_rate_22)
        call worker%store_model_field('Strain Rate-33', 'value', strain_rate_33)
        call worker%store_model_field('Strain Rate-12', 'value', strain_rate_12)
        call worker%store_model_field('Strain Rate-13', 'value', strain_rate_13)
        call worker%store_model_field('Strain Rate-23', 'value', strain_rate_23)
        call worker%store_model_field('Strain Rate-Trace', 'value', strain_rate_trace)
    end subroutine compute
    !***************************************************************************************




end module model_strain_rate
