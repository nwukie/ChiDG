module fluid_viscous_ale_volume_operator
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF

    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    implicit none

    private

    
    !> Volume flux for Fluid Viscous Terms.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/23/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: fluid_viscous_ale_volume_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type fluid_viscous_ale_volume_operator_t
    !******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/23/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(fluid_viscous_ale_volume_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name('Fluid Viscous ALE Volume Operator')

        ! Set operator type
        call self%set_operator_type('Volume Diffusive Operator')

        ! Set operator equations
        call self%add_primary_field('Density'   )
        call self%add_primary_field('Momentum-1')
        call self%add_primary_field('Momentum-2')
        call self%add_primary_field('Momentum-3')
        call self%add_primary_field('Energy'    )

    end subroutine init
    !********************************************************************************



    !> Volume flux routine for Fluid Viscous Terms.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/23/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(fluid_viscous_ale_volume_operator_t), intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:) ::    &
            density, mom1, mom2, mom3, energy,      &
            u, v, w, invdensity,                    &
            grad1_T, grad2_T, grad3_T,              &
            k,       k_l,     k_t,                  &
            tau_11, tau_22, tau_33,                 &
            tau_12, tau_13, tau_23,                 &
            flux_1, flux_2, flux_3

        real(rk),   allocatable, dimension(:)   :: r

        type(AD_D), allocatable                 :: flux_ref(:,:)

        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_primary_field_value_ale_element('Density'   )
        mom1    = worker%get_primary_field_value_ale_element('Momentum-1')
        mom2    = worker%get_primary_field_value_ale_element('Momentum-2')
        mom3    = worker%get_primary_field_value_ale_element('Momentum-3')
        energy  = worker%get_primary_field_value_ale_element('Energy'    )


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        r = worker%coordinate('1','volume')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2 = mom2 / r
        else if (worker%coordinate_system() == 'Cartesian') then

        else
            call chidg_signal(FATAL,"inlet, bad coordinate system")
        end if


        !
        ! Get Model fields:
        !   Second Coefficient of Viscosity
        !   Thermal Conductivity
        !
        k_l = worker%get_model_field_element('Laminar Thermal Conductivity',   'value')
        k_t = worker%get_model_field_element('Turbulent Thermal Conductivity', 'value')


        !
        ! compute effective conductivity. Laminar + Turbulent.
        !
        k = k_l + k_t


        !
        ! compute velocities
        !
        invdensity = ONE/density
        u          = mom1*invdensity
        v          = mom2*invdensity
        w          = mom3*invdensity


        !
        ! get temperature gradient
        !
        grad1_T = worker%get_model_field_element('Temperature Gradient - 1', 'value')
        grad2_T = worker%get_model_field_element('Temperature Gradient - 2', 'value')
        grad3_T = worker%get_model_field_element('Temperature Gradient - 3', 'value')


        !
        ! get shear stress components
        !
        tau_11 = worker%get_model_field_element('Shear-11', 'value')
        tau_22 = worker%get_model_field_element('Shear-22', 'value')
        tau_33 = worker%get_model_field_element('Shear-33', 'value')

        tau_12 = worker%get_model_field_element('Shear-12', 'value')
        tau_13 = worker%get_model_field_element('Shear-13', 'value')
        tau_23 = worker%get_model_field_element('Shear-23', 'value')





        !----------------------------------
        !            mass flux
        !----------------------------------


        !----------------------------------
        !         momentum-1 flux
        !----------------------------------
        flux_1 = -tau_11
        flux_2 = -tau_12
        flux_3 = -tau_13

        flux_ref = worker%post_process_volume_diffusive_flux_ale(flux_1, flux_2, flux_3)
        
        call worker%integrate_volume('Momentum-1',flux_ref(:,1),flux_ref(:,2),flux_ref(:,3))

        !----------------------------------
        !         momentum-2 flux
        !----------------------------------
        flux_1 = -tau_12
        flux_2 = -tau_22
        flux_3 = -tau_23

        flux_ref = worker%post_process_volume_diffusive_flux_ale(flux_1, flux_2, flux_3)
        !
        ! Convert to tangential to angular momentum flux
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            flux_1 = flux_1 * r
            flux_2 = flux_2 * r
            flux_3 = flux_3 * r
        else if (worker%coordinate_system() == 'Cartesian') then

        else
            call chidg_signal(FATAL,"inlet, bad coordinate system")
        end if

        call worker%integrate_volume('Momentum-2',flux_ref(:,1),flux_ref(:,2),flux_ref(:,3))

        !----------------------------------
        !         momentum-3 flux
        !----------------------------------
        flux_1 = -tau_13
        flux_2 = -tau_23
        flux_3 = -tau_33

        flux_ref = worker%post_process_volume_diffusive_flux_ale(flux_1, flux_2, flux_3)
        call worker%integrate_volume('Momentum-3',flux_ref(:,1),flux_ref(:,2),flux_ref(:,3))

        !----------------------------------
        !           energy flux
        !----------------------------------
        flux_1 = -k*grad1_T  -  (u*tau_11 + v*tau_12 + w*tau_13)
        flux_2 = -k*grad2_T  -  (u*tau_12 + v*tau_22 + w*tau_23)
        flux_3 = -k*grad3_T  -  (u*tau_13 + v*tau_23 + w*tau_33)

        flux_ref = worker%post_process_volume_diffusive_flux_ale(flux_1, flux_2, flux_3)
        call worker%integrate_volume('Energy',flux_ref(:,1),flux_ref(:,2),flux_ref(:,3))

    end subroutine compute
    !*********************************************************************************************************






end module fluid_viscous_ale_volume_operator
