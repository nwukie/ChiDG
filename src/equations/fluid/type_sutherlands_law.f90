module type_sutherlands_law
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: THREE, TWO
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !>  Sutherland's Law for Laminary Viscosity.
    !!
    !!  Model Fields:
    !!      - Viscosity
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/3/2016
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: sutherlands_law_t

    contains

        procedure   :: init
        procedure   :: compute

    end type sutherlands_law_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(sutherlands_law_t), intent(inout)   :: self

        call self%set_name('Sutherlands Law')

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
        class(sutherlands_law_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            viscosity, T

        real(rk) :: mu0 = 1.7894e-5_rk  ! [kg/(m*s)]
        real(rk) :: T0  = 273.11_rk     ! [K]
        real(rk) :: S   = 110.56_rk     ! [K]

        !
        ! Interpolate solution to quadrature nodes
        !
        T = worker%get_model_field_general('Temperature','value')


        !
        ! Sutherlands Law for Laminar Viscosity
        !
        viscosity = mu0*((T/T0)**(THREE/TWO))*(T0+S)/(T+S)


        !
        ! Contribute laminar viscosity
        !
        call worker%store_model_field('Laminar Viscosity', 'value', viscosity)


    end subroutine compute
    !***************************************************************************************




end module type_sutherlands_law
