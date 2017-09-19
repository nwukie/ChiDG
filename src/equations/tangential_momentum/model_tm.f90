module model_tm
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, ZERO, PI, TWO
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
    type, extends(model_t)  :: tm_t

    contains

        procedure   :: init
        procedure   :: compute

    end type tm_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(tm_t), intent(inout)   :: self

        call self%set_name('TM')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('Density')
        call self%add_model_field('Velocity-1')
        call self%add_model_field('Velocity-2')


    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Pressure and Temperature.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(tm_t),           intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            p, density, u, v

        real(rk),   dimension(:),   allocatable :: &
            r, theta

        real(rk) :: slope

        !
        ! Interpolate solution to quadrature nodes
        !
        p = worker%get_field('Pressure', 'value')

        !
        ! Initialize storage
        !
        density = p
        u       = p
        v       = p
        density = ZERO
        u       = ZERO
        v       = ZERO
        


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        r     = worker%coordinate('1')
        theta = worker%coordinate('2')


        slope   = 10._rk    ! (m/s) / m
        density = 1.19_rk
        u       = 0._rk
        !u       = 10._rk + r*slope
        !v       = 10._rk*sin(4._rk*theta)
        v       = 10._rk + r*slope
        !v = ZERO


        call worker%store_model_field('Density',    'value', density)
        call worker%store_model_field('Velocity-1', 'value', u      )
        call worker%store_model_field('Velocity-2', 'value', v      )


    end subroutine compute
    !***************************************************************************************









end module model_tm
