module fluid_viscous_volume_operator
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
    type, extends(operator_t), public :: fluid_viscous_volume_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type fluid_viscous_volume_operator_t
    !******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/23/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(fluid_viscous_volume_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name('Fluid Viscous Volume Operator')

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
        class(fluid_viscous_volume_operator_t), intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:) ::    &
            density, mom1, mom2, mom3, energy,      &
            u, v, w, invdensity,                    &
            grad1_T, grad2_T, grad3_T,              &
            k,       k_l,     k_t,                  &
            tau_11, tau_22, tau_33,                 &
            tau_12, tau_13, tau_23,                 &
            flux_1, flux_2, flux_3, r

        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_field('Density'   , 'value', 'element')
        mom1    = worker%get_field('Momentum-1', 'value', 'element')
        mom2    = worker%get_field('Momentum-2', 'value', 'element')
        mom3    = worker%get_field('Momentum-3', 'value', 'element')
        energy  = worker%get_field('Energy'    , 'value', 'element')


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','volume')
            mom2 = mom2 / r
        end if


        !
        ! Get Model fields:
        !   Second Coefficient of Viscosity
        !   Thermal Conductivity
        !
        k_l = worker%get_field('Laminar Thermal Conductivity',   'value', 'element')
        k_t = worker%get_field('Turbulent Thermal Conductivity', 'value', 'element')


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
        grad1_T = worker%get_field('Temperature Gradient - 1', 'value', 'element')
        grad2_T = worker%get_field('Temperature Gradient - 2', 'value', 'element')
        grad3_T = worker%get_field('Temperature Gradient - 3', 'value', 'element')


        !
        ! get shear stress components
        !
        tau_11 = worker%get_field('Shear-11', 'value', 'element')
        tau_22 = worker%get_field('Shear-22', 'value', 'element')
        tau_33 = worker%get_field('Shear-33', 'value', 'element')

        tau_12 = worker%get_field('Shear-12', 'value', 'element')
        tau_13 = worker%get_field('Shear-13', 'value', 'element')
        tau_23 = worker%get_field('Shear-23', 'value', 'element')





        !----------------------------------
        !            mass flux
        !----------------------------------


        !----------------------------------
        !         momentum-1 flux
        !----------------------------------
        flux_1 = -tau_11
        flux_2 = -tau_12
        flux_3 = -tau_13
        
        call worker%integrate_volume_flux('Momentum-1','Diffusion',flux_1,flux_2,flux_3)

        !----------------------------------
        !         momentum-2 flux
        !----------------------------------
        flux_1 = -tau_12
        flux_2 = -tau_22
        flux_3 = -tau_23

        ! Convert to tangential to angular momentum flux
        if (worker%coordinate_system() == 'Cylindrical') then
            flux_1 = flux_1 * r
            flux_2 = flux_2 * r
            flux_3 = flux_3 * r
        end if

        call worker%integrate_volume_flux('Momentum-2','Diffusion',flux_1,flux_2,flux_3)

        !----------------------------------
        !         momentum-3 flux
        !----------------------------------
        flux_1 = -tau_13
        flux_2 = -tau_23
        flux_3 = -tau_33

        call worker%integrate_volume_flux('Momentum-3','Diffusion',flux_1,flux_2,flux_3)

        !----------------------------------
        !           energy flux
        !----------------------------------
        flux_1 = -k*grad1_T  -  (u*tau_11 + v*tau_12 + w*tau_13)
        flux_2 = -k*grad2_T  -  (u*tau_12 + v*tau_22 + w*tau_23)
        flux_3 = -k*grad3_T  -  (u*tau_13 + v*tau_23 + w*tau_33)

        call worker%integrate_volume_flux('Energy','Diffusion',flux_1,flux_2,flux_3)

    end subroutine compute
    !***********************************************************************************






end module fluid_viscous_volume_operator
