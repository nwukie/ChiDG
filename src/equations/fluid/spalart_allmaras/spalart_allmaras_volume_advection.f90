module spalart_allmaras_volume_advection
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF
    use mod_fluid,              only: omega
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    type, extends(operator_t), public :: spalart_allmaras_volume_advection_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type spalart_allmaras_volume_advection_operator_t
    !***********************************************************************************************

contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine init(self)
        class(spalart_allmaras_volume_advection_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('Spalart-Allmaras Volume Advection Operator')

        ! Set operator type
        call self%set_operator_type('Volume Advective Operator')

        ! Set operator equations
        call self%add_primary_field('Density * NuTilde')

    end subroutine init
    !***********************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(spalart_allmaras_volume_advection_operator_t), intent(inout)   :: self
        type(chidg_worker_t),                                intent(inout)   :: worker
        class(properties_t),                                 intent(inout)   :: prop


        type(AD_D), dimension(:), allocatable   ::  &
            density_nutilde, u_a, v_a, w_a,         &
            flux_1, flux_2, flux_3


        !
        ! Interpolate solution to quadrature nodes
        !
        density_nutilde = worker%get_primary_field_element('Density * NuTilde', 'value')

        
        !
        ! Get fluid advection velocity
        ! 
        u_a = worker%get_model_field_element('Advection Velocity-1', 'value')
        v_a = worker%get_model_field_element('Advection Velocity-2', 'value')
        w_a = worker%get_model_field_element('Advection Velocity-3', 'value')


        !
        ! Compute average flux and field difference.
        ! 
        flux_1 = density_nutilde*u_a
        flux_2 = density_nutilde*v_a
        flux_3 = density_nutilde*w_a


        !
        ! Integrate flux
        !
        call worker%integrate_volume('Density * NuTilde',flux_1,flux_2,flux_3)


    end subroutine compute
    !************************************************************************************************









end module spalart_allmaras_volume_advection
