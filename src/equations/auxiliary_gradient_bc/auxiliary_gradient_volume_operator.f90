module auxiliary_gradient_volume_operator
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    use ieee_arithmetic
    implicit none
    private

    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/5/2018
    !!
    !-------------------------------------------------------------------------
    type, extends(operator_t), public :: auxiliary_gradient_volume_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type auxiliary_gradient_volume_operator_t
    !*************************************************************************

contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/5/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(auxiliary_gradient_volume_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name('Auxiliary Gradient Volume Operator')

        ! Set operator type
        call self%set_operator_type('Volume Diffusive Operator')

        ! Set operator equations
        call self%add_primary_field('Density_TEMP')
        call self%add_primary_field('Velocity-1_TEMP')
        call self%add_primary_field('Pressure_TEMP')

    end subroutine init
    !********************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/5/2018
    !!
    !------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(auxiliary_gradient_volume_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:)   ::                      &
            grad1_density_bc,   grad2_density_bc,   grad3_density_bc,   &
            grad1_vel1_bc,      grad2_vel1_bc,      grad3_vel1_bc,      &
            grad1_p_bc,         grad2_p_bc,         grad3_p_bc,         &
            grad1_density,      grad2_density,      grad3_density,      &
            grad1_vel1,         grad2_vel1,         grad3_vel1,         &
            grad1_p,            grad2_p,            grad3_p,            &
            flux_1, flux_2, flux_3, source
            

        !----------------------------------
        ! Auxiliary Velocity-1
        !----------------------------------
        grad1_density_bc = worker%get_field('Density_TEMP', 'grad1', 'element')
        grad2_density_bc = worker%get_field('Density_TEMP', 'grad2', 'element')
        grad3_density_bc = worker%get_field('Density_TEMP', 'grad3', 'element')

        grad1_density = worker%get_field('Density', 'grad1', 'element', override_lift=.true.)
        grad2_density = worker%get_field('Density', 'grad2', 'element', override_lift=.true.)
        grad3_density = worker%get_field('Density', 'grad3', 'element', override_lift=.true.)

        ! Compute volume flux at quadrature nodes
        flux_1 = grad1_density_bc - grad1_density
        flux_2 = grad2_density_bc - grad2_density
        flux_3 = grad3_density_bc - grad3_density

        ! Integrate volume flux
        call worker%integrate_volume_flux('Density_TEMP','Diffusion',flux_1,flux_2,flux_3)


        !----------------------------------
        ! Auxiliary Velocity-1
        !----------------------------------
        grad1_vel1_bc = worker%get_field('Velocity-1_TEMP', 'grad1', 'element')
        grad2_vel1_bc = worker%get_field('Velocity-1_TEMP', 'grad2', 'element')
        grad3_vel1_bc = worker%get_field('Velocity-1_TEMP', 'grad3', 'element')

        grad1_vel1 = worker%get_field('Velocity-1 Gradient-1', 'value', 'element')
        grad2_vel1 = worker%get_field('Velocity-1 Gradient-2', 'value', 'element')
        grad3_vel1 = worker%get_field('Velocity-1 Gradient-3', 'value', 'element')

        ! Compute volume flux at quadrature nodes
        flux_1 = grad1_vel1_bc - grad1_vel1
        flux_2 = grad2_vel1_bc - grad2_vel1
        flux_3 = grad3_vel1_bc - grad3_vel1

        ! Integrate volume flux
        call worker%integrate_volume_flux('Velocity-1_TEMP','Diffusion',flux_1,flux_2,flux_3)




        !----------------------------------
        ! Auxiliary Pressure
        !----------------------------------
        grad1_p_bc = worker%get_field('Pressure_TEMP', 'grad1', 'element')
        grad2_p_bc = worker%get_field('Pressure_TEMP', 'grad2', 'element')
        grad3_p_bc = worker%get_field('Pressure_TEMP', 'grad3', 'element')

        grad1_p = worker%get_field('Pressure Gradient - 1', 'value', 'element')
        grad2_p = worker%get_field('Pressure Gradient - 2', 'value', 'element')
        grad3_p = worker%get_field('Pressure Gradient - 3', 'value', 'element')

        ! Compute volume flux at quadrature nodes
        flux_1 = grad1_p_bc - grad1_p
        flux_2 = grad2_p_bc - grad2_p
        flux_3 = grad3_p_bc - grad3_p

        ! Integrate volume flux
        call worker%integrate_volume_flux('Pressure_TEMP','Diffusion',flux_1,flux_2,flux_3)



    end subroutine compute
    !****************************************************************************************************






end module auxiliary_gradient_volume_operator
