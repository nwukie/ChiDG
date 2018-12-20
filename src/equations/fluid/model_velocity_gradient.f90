module model_velocity_gradient
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, TWO
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    use ieee_arithmetic
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
    type, extends(model_t)  :: velocity_gradient_t


    contains

        procedure   :: init
        procedure   :: compute

    end type velocity_gradient_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(velocity_gradient_t), intent(inout)   :: self

        call self%set_name('Velocity Gradient')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Velocity 1 - Gradient 1')
        call self%add_model_field('Velocity 1 - Gradient 2')
        call self%add_model_field('Velocity 1 - Gradient 3')
        call self%add_model_field('Velocity 2 - Gradient 1')
        call self%add_model_field('Velocity 2 - Gradient 2')
        call self%add_model_field('Velocity 2 - Gradient 3')
        call self%add_model_field('Velocity 3 - Gradient 1')
        call self%add_model_field('Velocity 3 - Gradient 2')
        call self%add_model_field('Velocity 3 - Gradient 3')

    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Pressure and Temperature.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(velocity_gradient_t),  intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), allocatable,    dimension(:) ::         &
            density, mom1, mom2, mom3,                      &
            grad1_density, grad2_density, grad3_density,    &
            grad1_mom1,    grad2_mom1,    grad3_mom1,       &
            grad1_mom2,    grad2_mom2,    grad3_mom2,       &
            grad1_mom3,    grad2_mom3,    grad3_mom3,       &
            grad1_u,       grad2_u,       grad3_u,          &
            grad1_v,       grad2_v,       grad3_v,          &
            grad1_w,       grad2_w,       grad3_w,          &
            du_ddensity,   dv_ddensity,   dw_ddensity,      &
            du_dmom1,      dv_dmom2,      dw_dmom3,         &
            invdensity

        real(rk),   allocatable,    dimension(:) :: r


        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_field('Density',    'value')
        mom1    = worker%get_field('Momentum-1', 'value')
        mom2    = worker%get_field('Momentum-2', 'value')
        mom3    = worker%get_field('Momentum-3', 'value')

        if (any(ieee_is_nan(density(:)%x_ad_))) print *, 'velgrad: density is nan'
        if (any(ieee_is_nan(mom1(:)%x_ad_))) print *, 'velgrad: mom1 is nan'
        if (any(ieee_is_nan(mom2(:)%x_ad_))) print *, 'velgrad: mom2 is nan'
        if (any(ieee_is_nan(mom3(:)%x_ad_))) print *, 'velgrad: mom3 is nan'

        !
        ! Interpolate gradient to quadrature nodes
        !
        grad1_density    = worker%get_field('Density'   , 'grad1')!, user_no_lift=.true.)
        grad2_density    = worker%get_field('Density'   , 'grad2')!, user_no_lift=.true.)
        grad3_density    = worker%get_field('Density'   , 'grad3')!, user_no_lift=.true.)

        grad1_mom1       = worker%get_field('Momentum-1', 'grad1')!, user_no_lift=.true.)
        grad2_mom1       = worker%get_field('Momentum-1', 'grad2')!, user_no_lift=.true.)
        grad3_mom1       = worker%get_field('Momentum-1', 'grad3')!, user_no_lift=.true.)

        grad1_mom2       = worker%get_field('Momentum-2', 'grad1')!, user_no_lift=.true.)
        grad2_mom2       = worker%get_field('Momentum-2', 'grad2')!, user_no_lift=.true.)
        grad3_mom2       = worker%get_field('Momentum-2', 'grad3')!, user_no_lift=.true.)

        grad1_mom3       = worker%get_field('Momentum-3', 'grad1')!, user_no_lift=.true.)
        grad2_mom3       = worker%get_field('Momentum-3', 'grad2')!, user_no_lift=.true.)
        grad3_mom3       = worker%get_field('Momentum-3', 'grad3')!, user_no_lift=.true.)

        if (any(ieee_is_nan(grad1_density(:)%x_ad_))) print *, 'velgrad: grad1_density is nan'
        if (any(ieee_is_nan(grad2_density(:)%x_ad_))) print *, 'velgrad: grad2_density is nan'
        if (any(ieee_is_nan(grad3_density(:)%x_ad_))) print *, 'velgrad: grad3_density is nan'


        if (any(ieee_is_nan(grad1_mom1(:)%x_ad_))) print *, 'velgrad: grad1_mom1 is nan'
        if (any(ieee_is_nan(grad2_mom1(:)%x_ad_))) print *, 'velgrad: grad2_mom1 is nan'
        if (any(ieee_is_nan(grad3_mom1(:)%x_ad_))) print *, 'velgrad: grad3_mom1 is nan'


        if (any(ieee_is_nan(grad1_mom2(:)%x_ad_))) print *, 'velgrad: grad1_mom2 is nan'
        if (any(ieee_is_nan(grad2_mom2(:)%x_ad_))) print *, 'velgrad: grad2_mom2 is nan'
        if (any(ieee_is_nan(grad3_mom2(:)%x_ad_))) print *, 'velgrad: grad3_mom2 is nan'

        if (any(ieee_is_nan(grad1_mom3(:)%x_ad_))) print *, 'velgrad: grad1_mom3 is nan'
        if (any(ieee_is_nan(grad2_mom3(:)%x_ad_))) print *, 'velgrad: grad2_mom3 is nan'
        if (any(ieee_is_nan(grad3_mom3(:)%x_ad_))) print *, 'velgrad: grad3_mom3 is nan'


        
        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        ! Also convert derivatives from derivatives of angular momentum to tangential.
        !
        ! We want:
        !       (rho * u_theta)  instead of      (r * rho * u_theta)
        !   grad(rho * u_theta)  instead of  grad(r * rho * u_theta)
        !
        !   grad(rho * u_theta) = grad(r * rho * u_theta)/r  -  grad(r)(rho*u_theta)/r
        !
        ! Where grad(r) = [1,0,0]
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1')
            mom2       = mom2 / r
            grad1_mom2 = (grad1_mom2/r) - mom2/r
            grad2_mom2 = (grad2_mom2/r)
            grad3_mom2 = (grad3_mom2/r)
        end if



        invdensity = ONE/(density )
        if (any(ieee_is_nan(invdensity(:)%x_ad_))) print*, 'inv density is nan in velgrad'


        
        !
        ! compute velocity jacobians
        !
        du_ddensity  = -invdensity*invdensity*mom1
        dv_ddensity  = -invdensity*invdensity*mom2
        dw_ddensity  = -invdensity*invdensity*mom3

        du_dmom1 = invdensity
        dv_dmom2 = invdensity
        dw_dmom3 = invdensity



        !
        ! compute velocity gradients via chain rule:
        !
        !   u = f(rho,rhou)
        !
        !   grad(u) = dudrho * grad(rho)  +  dudrhou * grad(rhou)
        !
        grad1_u = du_ddensity*grad1_density  +  du_dmom1*grad1_mom1
        grad2_u = du_ddensity*grad2_density  +  du_dmom1*grad2_mom1
        grad3_u = du_ddensity*grad3_density  +  du_dmom1*grad3_mom1

        grad1_v = dv_ddensity*grad1_density  +  dv_dmom2*grad1_mom2
        grad2_v = dv_ddensity*grad2_density  +  dv_dmom2*grad2_mom2
        grad3_v = dv_ddensity*grad3_density  +  dv_dmom2*grad3_mom2

        grad1_w = dw_ddensity*grad1_density  +  dw_dmom3*grad1_mom3
        grad2_w = dw_ddensity*grad2_density  +  dw_dmom3*grad2_mom3
        grad3_w = dw_ddensity*grad3_density  +  dw_dmom3*grad3_mom3

        if (any(ieee_is_nan(grad1_u(:)%x_ad_))) print *, 'grad1_u is nan'
        if (any(ieee_is_nan(grad2_u(:)%x_ad_))) print *, 'grad2_u is nan'
        if (any(ieee_is_nan(grad3_u(:)%x_ad_))) print *, 'grad3_u is nan'


        if (any(ieee_is_nan(grad1_v(:)%x_ad_))) print *, 'grad1_v is nan'
        if (any(ieee_is_nan(grad2_v(:)%x_ad_))) print *, 'grad2_v is nan'
        if (any(ieee_is_nan(grad3_v(:)%x_ad_))) print *, 'grad3_v is nan'

        if (any(ieee_is_nan(grad1_w(:)%x_ad_))) print *, 'grad1_w is nan'
        if (any(ieee_is_nan(grad2_w(:)%x_ad_))) print *, 'grad2_w is nan'
        if (any(ieee_is_nan(grad3_w(:)%x_ad_))) print *, 'grad3_w is nan'




        call worker%store_model_field('Velocity 1 - Gradient 1', 'value', grad1_u)
        call worker%store_model_field('Velocity 1 - Gradient 2', 'value', grad2_u)
        call worker%store_model_field('Velocity 1 - Gradient 3', 'value', grad3_u)

        call worker%store_model_field('Velocity 2 - Gradient 1', 'value', grad1_v)
        call worker%store_model_field('Velocity 2 - Gradient 2', 'value', grad2_v)
        call worker%store_model_field('Velocity 2 - Gradient 3', 'value', grad3_v)

        call worker%store_model_field('Velocity 3 - Gradient 1', 'value', grad1_w)
        call worker%store_model_field('Velocity 3 - Gradient 2', 'value', grad2_w)
        call worker%store_model_field('Velocity 3 - Gradient 3', 'value', grad3_w)

    end subroutine compute
    !***************************************************************************************




end module model_velocity_gradient
