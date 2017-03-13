module type_temperature_gradient
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
    !!      : shear_21, shear_22, shear_23
    !!      : shear_31, shear_32, shear_33
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: temperature_gradient_t

        real(rk)    :: gam = 1.4_rk     ! ratio of specific heats
        real(rk)    :: R   = 287.15_rk  ! ideal gas constant [J/(kg*K)]

    contains

        procedure   :: init
        procedure   :: compute

    end type temperature_gradient_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(temperature_gradient_t), intent(inout)   :: self

        call self%set_name('Temperature Gradient')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Temperature Gradient - 1')
        call self%add_model_field('Temperature Gradient - 2')
        call self%add_model_field('Temperature Gradient - 3')

    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Pressure and Temperature.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(temperature_gradient_t),  intent(in)      :: self
        type(chidg_worker_t),           intent(inout)   :: worker

        type(AD_D), allocatable,    dimension(:) ::         &
            density, mom1, mom2, mom3, energy, invdensity,  &
            grad1_density, grad2_density, grad3_density,    &
            grad1_mom1,    grad2_mom1,    grad3_mom1,       &
            grad1_mom2,    grad2_mom2,    grad3_mom2,       &
            grad1_mom3,    grad2_mom3,    grad3_mom3,       &
            grad1_energy,  grad2_energy,  grad3_energy,     &
            grad1_T,       grad2_T,       grad3_T,          &
            u,             v,             w,       p,       &
            dT_ddensity,   dT_dmom1,   dT_dmom2,   dT_dmom3,   dT_denergy,  &
            dp_ddensity,   dp_dmom1,   dp_dmom2,   dp_dmom3,   dp_denergy,  &
            dke_ddensity,  dke_dmom1,  dke_dmom2,  dke_dmom3




        real(rk),   allocatable,    dimension(:) :: r
        real(rk)    :: gam, const


        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_primary_field_general('Density',    'value')
        mom1    = worker%get_primary_field_general('Momentum-1', 'value')
        mom2    = worker%get_primary_field_general('Momentum-2', 'value')
        mom3    = worker%get_primary_field_general('Momentum-3', 'value')
        energy  = worker%get_primary_field_general('Energy',     'value')


        !
        ! Interpolate gradient to quadrature nodes
        !
        grad1_density = worker%get_primary_field_general('Density'   , 'grad1+lift')
        grad2_density = worker%get_primary_field_general('Density'   , 'grad2+lift')
        grad3_density = worker%get_primary_field_general('Density'   , 'grad3+lift')

        grad1_mom1    = worker%get_primary_field_general('Momentum-1', 'grad1+lift')
        grad2_mom1    = worker%get_primary_field_general('Momentum-1', 'grad2+lift')
        grad3_mom1    = worker%get_primary_field_general('Momentum-1', 'grad3+lift')

        grad1_mom2    = worker%get_primary_field_general('Momentum-2', 'grad1+lift')
        grad2_mom2    = worker%get_primary_field_general('Momentum-2', 'grad2+lift')
        grad3_mom2    = worker%get_primary_field_general('Momentum-2', 'grad3+lift')

        grad1_mom3    = worker%get_primary_field_general('Momentum-3', 'grad1+lift')
        grad2_mom3    = worker%get_primary_field_general('Momentum-3', 'grad2+lift')
        grad3_mom3    = worker%get_primary_field_general('Momentum-3', 'grad3+lift')

        grad1_energy  = worker%get_primary_field_general('Energy'    , 'grad1+lift')
        grad2_energy  = worker%get_primary_field_general('Energy'    , 'grad2+lift')
        grad3_energy  = worker%get_primary_field_general('Energy'    , 'grad3+lift')



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












        !
        ! Get model field 'Pressure'
        !
        p = worker%get_model_field_general('Pressure', 'value')
        gam = 1.4_rk




        !
        ! Compute velocities
        !
        invdensity = ONE/density
        u          = mom1*invdensity
        v          = mom2*invdensity
        w          = mom3*invdensity


        !
        ! Compute Kinetic Energy Jacobians
        !
        dke_ddensity = -HALF*(u*u + v*v + w*w)
        dke_dmom1    = u
        dke_dmom2    = v
        dke_dmom3    = w



        !
        ! Compute Pressure Jacobians
        !
        dp_ddensity = -(gam-ONE)*dke_ddensity
        dp_dmom1    = -(gam-ONE)*dke_dmom1
        dp_dmom2    = -(gam-ONE)*dke_dmom2
        dp_dmom3    = -(gam-ONE)*dke_dmom3
        dp_denergy  =  dp_dmom3    ! Initialize derivatives
        dp_denergy  =  (gam-ONE)   ! No negative sign


        !
        ! Compute Temperature Jacobians
        !
        const = ONE/287.15_rk
        dT_ddensity = const*invdensity*dp_ddensity  -  const*invdensity*invdensity*p
        dT_dmom1    = const*invdensity*dp_dmom1
        dT_dmom2    = const*invdensity*dp_dmom2
        dT_dmom3    = const*invdensity*dp_dmom3
        dT_denergy  = const*invdensity*dp_denergy


        !
        ! Compute temperature gradient
        !
        grad1_T = dT_ddensity*grad1_density + dT_dmom1*grad1_mom1 + dT_dmom2*grad1_mom2 + dT_dmom3*grad1_mom3 + dT_denergy*grad1_energy
        grad2_T = dT_ddensity*grad2_density + dT_dmom1*grad2_mom1 + dT_dmom2*grad2_mom2 + dT_dmom3*grad2_mom3 + dT_denergy*grad2_energy
        grad3_T = dT_ddensity*grad3_density + dT_dmom1*grad3_mom1 + dT_dmom2*grad3_mom2 + dT_dmom3*grad3_mom3 + dT_denergy*grad3_energy



        call worker%store_model_field('Temperature Gradient - 1', 'value', grad1_T)
        call worker%store_model_field('Temperature Gradient - 2', 'value', grad2_T)
        call worker%store_model_field('Temperature Gradient - 3', 'value', grad3_T)


    end subroutine compute
    !***************************************************************************************




end module type_temperature_gradient
