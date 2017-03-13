module type_reynolds_analogy
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: THREE, TWO
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !>  Reynolds' Analogy used to compute thermal conductivity.
    !!
    !!  Model Fields:
    !!      - Thermal Conductivity
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: reynolds_analogy_t

        real(rk)    :: Cp = 1003.0_rk
        real(rk)    :: Pr = 0.8_rk

    contains

        procedure   :: init
        procedure   :: compute

    end type reynolds_analogy_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(reynolds_analogy_t), intent(inout)   :: self

        call self%set_name('Reynolds Analogy')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('Laminar Thermal Conductivity')


    end subroutine init
    !***************************************************************************************






    !>  Routine for computing a viscosity contribution from Sutherland's Law.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(reynolds_analogy_t),  intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: viscosity, thermal_conductivity


        !
        ! Interpolate solution to quadrature nodes
        !
        viscosity = worker%get_model_field_general('Laminar Viscosity','value')


        !
        ! Stokes' Hypothesis for the second coefficient of viscosity
        !
        thermal_conductivity = self%Cp * viscosity / self%Pr


        !
        ! Contribute second coefficient of viscosity
        !
        call worker%store_model_field('Laminar Thermal Conductivity', 'value', thermal_conductivity)


    end subroutine compute
    !***************************************************************************************




end module type_reynolds_analogy
