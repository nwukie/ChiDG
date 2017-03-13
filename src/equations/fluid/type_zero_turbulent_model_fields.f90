module type_zero_turbulent_model_fields
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: ZERO
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    implicit none


    


    !>  Provide zero-values for turbulent model fields, implying laminar flow.
    !!
    !!  Model Fields:
    !!      - Turbulent Viscosity
    !!      - Second Coefficient of Turbulent Viscosity
    !!      - Turbulent Thermal Conductivity
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: zero_turbulent_model_fields_t

    contains

        procedure   :: init
        procedure   :: compute

    end type zero_turbulent_model_fields_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(zero_turbulent_model_fields_t), intent(inout)   :: self

        call self%set_name('Zero Turbulent Model Fields')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('Turbulent Viscosity')
        call self%add_model_field('Second Coefficient of Turbulent Viscosity')
        call self%add_model_field('Turbulent Thermal Conductivity')


    end subroutine init
    !***************************************************************************************






    !>  Routine for computing a viscosity contribution from Sutherland's Law.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/1/2016
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(zero_turbulent_model_fields_t), intent(in)      :: self
        type(chidg_worker_t),                   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: vals


        !
        ! Get field to set derivative arrays. Set mut to zero.
        !
        vals = worker%get_primary_field_general('Density','value')
        vals = ZERO

        !
        ! Contribute second coefficient of viscosity
        !
        call worker%store_model_field('Turbulent Viscosity',                       'value', vals)
        call worker%store_model_field('Second Coefficient of Turbulent Viscosity', 'value', vals)
        call worker%store_model_field('Turbulent Thermal Conductivity',            'value', vals)


    end subroutine compute
    !***************************************************************************************




end module type_zero_turbulent_model_fields
