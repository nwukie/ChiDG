module model_velocity_div_curl
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
    type, extends(model_t)  :: velocity_div_curl_t


    contains

        procedure   :: init
        procedure   :: compute

    end type velocity_div_curl_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(velocity_div_curl_t), intent(inout)   :: self

        call self%set_name('Velocity Divergence and Curl')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Velocity Divergence')
        call self%add_model_field('Velocity Curl - 1')
        call self%add_model_field('Velocity Curl - 2')
        call self%add_model_field('Velocity Curl - 3')
    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Pressure and Temperature.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(velocity_div_curl_t),  intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), allocatable,    dimension(:) ::         &
            grad1_u,       grad2_u,       grad3_u,          &
            grad1_v,       grad2_v,       grad3_v,          &
            grad1_w,       grad2_w,       grad3_w,          &
            curl_vel_1, curl_vel_2, curl_vel_3,             &
            div_vel 

        real(rk),   allocatable,    dimension(:) :: r



         grad1_u = worker%get_field('Velocity 1 - Gradient 1', 'value')
         grad2_u = worker%get_field('Velocity 1 - Gradient 2', 'value')
         grad3_u = worker%get_field('Velocity 1 - Gradient 3', 'value')
                  
         grad1_v = worker%get_field('Velocity 2 - Gradient 1', 'value')
         grad2_v = worker%get_field('Velocity 2 - Gradient 2', 'value')
         grad3_v = worker%get_field('Velocity 2 - Gradient 3', 'value')
                  
         grad1_w = worker%get_field('Velocity 3 - Gradient 1', 'value')
         grad2_w = worker%get_field('Velocity 3 - Gradient 2', 'value')
         grad3_w = worker%get_field('Velocity 3 - Gradient 3', 'value')

         div_vel = grad1_u + grad2_v + grad3_w

         curl_vel_1 = grad2_w-grad3_v
         curl_vel_2 = grad3_u-grad1_w
         curl_vel_3 = grad1_v-grad2_u

         call worker%store_model_field('Velocity Divergence', 'value', div_vel)

         call worker%store_model_field('Velocity Curl - 1', 'value', curl_vel_1)
         call worker%store_model_field('Velocity Curl - 2', 'value', curl_vel_2)
         call worker%store_model_field('Velocity Curl - 3', 'value', curl_vel_3)

    end subroutine compute
    !***************************************************************************************




end module model_velocity_div_curl
