module model_entropy
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, ZERO
    use mod_fluid,          only: Rgas, gam
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
    type, extends(model_t)  :: entropy_t

    contains

        procedure   :: init
        procedure   :: compute

    end type entropy_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(entropy_t), intent(inout)   :: self

        call self%set_name('Entropy')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('Entropy')


    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Pressure and Temperature.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(entropy_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            density, pressure, entropy 


        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_field('Density',    'value')
        pressure  = worker%get_field('Pressure',     'value')



        entropy = pressure/(density**gam)



        call worker%store_model_field('Entropy',    'value', entropy)



    end subroutine compute
    !***************************************************************************************









end module model_entropy
