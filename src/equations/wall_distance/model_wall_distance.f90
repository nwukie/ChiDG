module model_wall_distance
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, TWO
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
!    use eqn_wall_distance,  only: get_p_poisson_parameter
    use DNAD_D
    implicit none


    


    !>  An equation of state model for an ideal gas.
    !!
    !!  Model Fields:
    !!      - Wall Distance
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/5/2016
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: wall_distance_m

    contains

        procedure   :: init
        procedure   :: compute

    end type wall_distance_m
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/5/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(wall_distance_m), intent(inout)   :: self

        call self%set_name('Wall Distance : p-Poisson Normalization')

        call self%add_model_field('Wall Distance')

    end subroutine init
    !***************************************************************************************






    !>  Routine for computing the Wall Distance using a normalization of the output
    !!  from a p-Poisson equation based method.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/5/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(wall_distance_m),     intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            d, dddx, dddy, dddz, d_normalization, mag2, rho

        real(rk) :: p


        ! Get primary field to initialize derivatives
        rho  = worker%get_primary_field_general('Density', 'value')
        d    = rho
        dddx = rho
        dddy = rho
        dddz = rho


        !
        ! Interpolate solution to quadrature nodes
        !
        d    = worker%get_auxiliary_field_general('Wall Distance : p-Poisson', 'value')
        dddx = worker%get_auxiliary_field_general('Wall Distance : p-Poisson', 'grad1')
        dddy = worker%get_auxiliary_field_general('Wall Distance : p-Poisson', 'grad2')
        dddz = worker%get_auxiliary_field_general('Wall Distance : p-Poisson', 'grad3')


        !
        ! Compute wall distance normalization
        !
        !p = get_p_poisson_parameter()
        p = 6._rk
        mag2 = dddx*dddx + dddy*dddy + dddz*dddz
        d_normalization = (((p/(p-ONE))*d) + mag2**(p/TWO))**((p-ONE)/p) - mag2**((p-ONE)/TWO)


        !
        ! Store Wall Distance to model field
        !
        call worker%store_model_field('Wall Distance','value', d_normalization)


    end subroutine compute
    !***************************************************************************************




end module model_wall_distance
