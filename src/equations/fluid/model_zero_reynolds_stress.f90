module model_zero_reynolds_stress
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, TWO, ZERO
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !>  A model for computing reynolds stress
    !!
    !!  Model Fields:
    !!      : reynolds_11, reynolds_12, reynolds_13
    !!      :           reynolds_22, reynolds_23
    !!      :                     reynolds_33
    !!
    !!  Lower-triangular components are not computed because the tensor is symmetric.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: zero_reynolds_stress_t


    contains

        procedure   :: init
        procedure   :: compute

    end type zero_reynolds_stress_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(zero_reynolds_stress_t), intent(inout)   :: self

        call self%set_name('Zero Reynolds Stress')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('Reynolds-11')
        call self%add_model_field('Reynolds-22')
        call self%add_model_field('Reynolds-33')
        call self%add_model_field('Reynolds-12')
        call self%add_model_field('Reynolds-13')
        call self%add_model_field('Reynolds-23')

    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Pressure and Temperature.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(zero_reynolds_stress_t),  intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), allocatable,    dimension(:) ::         &
            density, mom1, mom2, mom3, energy,              &
            grad1_density, grad2_density, grad3_density,    &
            grad1_mom1,    grad2_mom1,    grad3_mom1,       &
            grad1_mom2,    grad2_mom2,    grad3_mom2,       &
            grad1_mom3,    grad2_mom3,    grad3_mom3,       &
            grad1_energy,  grad2_energy,  grad3_energy,     &
            grad1_u,       grad2_u,       grad3_u,          &
            grad1_v,       grad2_v,       grad3_v,          &
            grad1_w,       grad2_w,       grad3_w,          &
            mu,            mu_l,          mu_t,             &
            lamda,         lamda_l,       lamda_t,          &
            reynolds_11,      reynolds_22,      reynolds_33,         &
            reynolds_12,      reynolds_13,      reynolds_23,         &
            du_ddensity,   dv_ddensity,   dw_ddensity,      &
            du_dmom1,      dv_dmom2,      dw_dmom3,         &
            invdensity, div_velocity, u, v

        real(rk),   allocatable,    dimension(:) :: r


        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_field('Density',    'value')
        
        density = ZERO

        call worker%store_model_field('Reynolds-11', 'value', density)
        call worker%store_model_field('Reynolds-22', 'value', density)
        call worker%store_model_field('Reynolds-33', 'value', density)
        call worker%store_model_field('Reynolds-12', 'value', density)
        call worker%store_model_field('Reynolds-13', 'value', density)
        call worker%store_model_field('Reynolds-23', 'value', density)


    end subroutine compute
    !***************************************************************************************




end module model_zero_reynolds_stress
