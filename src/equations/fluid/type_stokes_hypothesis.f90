module type_stokes_hypothesis
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: THREE, TWO
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !>  Stokes' Hypothesis for computing the second coefficient of viscosity.
    !!
    !!  Model Fields:
    !!      - Second Coefficient of Viscosity
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/3/2016
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: stokes_hypothesis_t

    contains

        procedure   :: init
        procedure   :: compute

    end type stokes_hypothesis_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(stokes_hypothesis_t), intent(inout)   :: self

        call self%set_name('Stokes Hypothesis')

        call self%add_model_field('Second Coefficient of Laminar Viscosity')


    end subroutine init
    !***************************************************************************************






    !>  Routine for computing a viscosity contribution from Sutherland's Law.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(stokes_hypothesis_t), intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: viscosity, second_viscosity


        !
        ! Interpolate solution to quadrature nodes
        !
        viscosity = worker%get_model_field_general('Laminar Viscosity','value')


        !
        ! Stokes' Hypothesis for the second coefficient of viscosity
        !
        second_viscosity = -(TWO/THREE)*viscosity


        !
        ! Contribute second coefficient of viscosity
        !
        call worker%store_model_field('Second Coefficient of Laminar Viscosity', 'value', second_viscosity)


    end subroutine compute
    !***************************************************************************************




end module type_stokes_hypothesis
