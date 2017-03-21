module type_fluid_advection_velocity
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE
    use mod_fluid,          only: omega
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !>  An equation of state model for an ideal gas.
    !!
    !!  Model Fields:
    !!      - Pressure
    !!      - Temperature
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: fluid_advection_velocity_t

        real(rk)    :: gam = 1.4_rk     ! ratio of specific heats
        real(rk)    :: R   = 287.15_rk  ! ideal gas constant [J/(kg*K)]

    contains

        procedure   :: init
        procedure   :: compute

    end type fluid_advection_velocity_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(fluid_advection_velocity_t), intent(inout)   :: self

        call self%set_name('Fluid Advection Velocity')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('Advection Velocity-1')
        call self%add_model_field('Advection Velocity-2')
        call self%add_model_field('Advection Velocity-3')


    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Pressure and Temperature.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(fluid_advection_velocity_t),    intent(in)      :: self
        type(chidg_worker_t),           intent(inout)   :: worker

        type(AD_D), allocatable,    dimension(:)    ::  &
            density, mom1, mom2, mom3,                  &
            u,   v,   w,                                &
            u_a, v_a, w_a, invdensity

        real(rk),   allocatable,    dimension(:)    :: r

        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_primary_field_general('Density',    'value')
        mom1    = worker%get_primary_field_general('Momentum-1', 'value')
        mom2    = worker%get_primary_field_general('Momentum-2', 'value')
        mom3    = worker%get_primary_field_general('Momentum-3', 'value')


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        r = worker%coordinate('1')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2 = mom2 / r
        end if



        !
        ! Compute fluid velocities
        !
        invdensity = ONE/density
        u = mom1*invdensity
        v = mom2*invdensity
        w = mom3*invdensity
        


        !
        ! Compute advection velocities
        !
        u_a = u
        v_a = v - omega*r
        w_a = w



        call worker%store_model_field('Advection Velocity-1', 'value', u_a)
        call worker%store_model_field('Advection Velocity-2', 'value', v_a)
        call worker%store_model_field('Advection Velocity-3', 'value', w_a)


    end subroutine compute
    !***************************************************************************************









end module type_fluid_advection_velocity
