module model_rotation_rate
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
    type, extends(model_t)  :: rotation_rate_t


    contains

        procedure   :: init
        procedure   :: compute

    end type rotation_rate_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(rotation_rate_t), intent(inout)   :: self

        call self%set_name('Rotation Rate')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Rotation Rate-12')
        call self%add_model_field('Rotation Rate-13')
        call self%add_model_field('Rotation Rate-23')

    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Pressure and Temperature.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(rotation_rate_t),  intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), allocatable,    dimension(:) ::         &
            density, mom1, mom2, mom3,                      &
            grad1_density, grad2_density, grad3_density,    &
            grad1_vel1,    grad2_vel1,    grad3_vel1,       &
            grad1_vel2,    grad2_vel2,    grad3_vel2,       &
            grad1_vel3,    grad2_vel3,    grad3_vel3,       &
            rotation_rate_12, rotation_rate_13, rotation_rate_23

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


        rotation_rate_12 = HALF*(grad2_vel1 - grad1_vel2)
        rotation_rate_13 = HALF*(grad3_vel1 - grad1_vel3)
        rotation_rate_23 = HALF*(grad3_vel2 - grad2_vel3)


        call worker%store_model_field('Rotation Rate-12', 'value', rotation_rate_12)
        call worker%store_model_field('Rotation Rate-13', 'value', rotation_rate_13)
        call worker%store_model_field('Rotation Rate-23', 'value', rotation_rate_23)
    end subroutine compute
    !***************************************************************************************




end module model_rotation_rate
