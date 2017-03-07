module type_constant_viscosity_rans
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: THREE, TWO
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !>  Constant Viscosity for Laminary Viscosity.
    !!
    !!  Model Fields:
    !!      - Viscosity
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/26/2017
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: constant_viscosity_rans_t

    contains

        procedure   :: init
        procedure   :: compute

    end type constant_viscosity_rans_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   01/26/2017
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(constant_viscosity_rans_t), intent(inout)   :: self

        call self%set_name('Constant Viscosity RANS')
        call self%set_dependency('Q-')

        call self%add_model_field('Laminar Viscosity')

    end subroutine init
    !***************************************************************************************






    !>  Routine for computing a viscosity contribution from Sutherland's Law.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(constant_viscosity_rans_t),    intent(in)      :: self
        type(chidg_worker_t),           intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            viscosity, T

        real(rk) :: mu0 = 0.00018_rk  ! [kg/(m*s)]

        !
        ! Interpolate solution to quadrature nodes
        !
        T = worker%get_model_field_general('Temperature','value')
    

        !
        ! Constant Viscosity for Laminar Viscosity
        !   - initialize derivatives first...
        !
        viscosity = T
        viscosity = mu0


        !
        ! Contribute laminar viscosity
        !
        call worker%store_model_field('Laminar Viscosity', 'value', viscosity)


    end subroutine compute
    !***************************************************************************************




end module type_constant_viscosity_rans
