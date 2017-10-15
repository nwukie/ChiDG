module model_rac
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
    type, extends(model_t)  :: rac_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rac_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(rac_t), intent(inout)   :: self

        call self%set_name('RAC')
        call self%set_dependency('f(Q-)')

!        call self%add_model_field('Density')
!        call self%add_model_field('Velocity-1')
!        call self%add_model_field('Velocity-2')
!        call self%add_model_field('Velocity-1 : Grad1')
!        call self%add_model_field('Velocity-1 : Grad2')
!        call self%add_model_field('Velocity-2 : Grad1')
!        call self%add_model_field('Velocity-2 : Grad2')
!        call self%add_model_field('Momentum-1 : Grad1')
!        call self%add_model_field('Momentum-1 : Grad2')
!        call self%add_model_field('Momentum-2 : Grad1')
!        call self%add_model_field('Momentum-2 : Grad2')



    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Pressure and Temperature.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(rac_t),           intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            p, density, u, v,   &
            vel1_grad1, vel1_grad2, &
            vel2_grad1, vel2_grad2, &
            mom1_grad1, mom1_grad2, &
            mom2_grad1, mom2_grad2

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
        vel1_grad1 = p
        vel1_grad2 = p
        vel2_grad1 = p
        vel2_grad2 = p
        mom1_grad1 = p
        mom1_grad2 = p
        mom2_grad1 = p
        mom2_grad2 = p

        density = ZERO
        u       = ZERO
        v       = ZERO
        vel1_grad1 = ZERO
        vel1_grad2 = ZERO
        vel2_grad1 = ZERO
        vel2_grad2 = ZERO
        mom1_grad1 = ZERO
        mom1_grad2 = ZERO
        mom2_grad1 = ZERO
        mom2_grad2 = ZERO
        

        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        r     = worker%coordinate('1')
        theta = worker%coordinate('2')


        slope   = 10._rk    ! (m/s) / m
        density = 1._rk
        u       = 0._rk

        !v          = 10._rk + r*slope
        !vel2_grad1 = 10._rk
        !mom2_grad1 = density*10._rk


        v          = 10._rk*sin(4._rk*theta)
        vel2_grad2 = (ONE/r)*(4._rk*10._rk)*cos(4._rk*theta)
        mom2_grad2 = (ONE/r)*(4._rk*10._rk)*cos(4._rk*theta)

        !v          = (ONE/r)*10._rk*sin(4._rk*theta)
        !vel2_grad1 = log(r)*10._rk*sin(4._rk*theta)
        !mom2_grad1 = log(r)*10._rk*sin(4._rk*theta)
        !vel2_grad2 = (ONE/(r*r))*(4._rk*10._rk)*cos(4._rk*theta)
        !mom2_grad2 = (ONE/(r*r))*(4._rk*10._rk)*cos(4._rk*theta)

        call worker%store_model_field('Density',    'value', density)
        call worker%store_model_field('Velocity-1', 'value', u      )
        call worker%store_model_field('Velocity-2', 'value', v      )

        call worker%store_model_field('Velocity-1 : Grad1', 'value', vel1_grad1)
        call worker%store_model_field('Velocity-1 : Grad2', 'value', vel1_grad2)
        call worker%store_model_field('Velocity-2 : Grad1', 'value', vel2_grad1)
        call worker%store_model_field('Velocity-2 : Grad2', 'value', vel2_grad2)
        call worker%store_model_field('Momentum-1 : Grad1', 'value', mom1_grad1)
        call worker%store_model_field('Momentum-1 : Grad2', 'value', mom1_grad2)
        call worker%store_model_field('Momentum-2 : Grad1', 'value', mom2_grad1)
        call worker%store_model_field('Momentum-2 : Grad2', 'value', mom2_grad2)


    end subroutine compute
    !***************************************************************************************









end module model_rac
