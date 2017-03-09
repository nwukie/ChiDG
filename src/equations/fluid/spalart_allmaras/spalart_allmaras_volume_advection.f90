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


        type(AD_D), dimension(:), allocatable   ::      &
            density, mom1, mom2, mom3, density_nutilde, &
            invdensity, u, v, w, flux_1, flux_2, flux_3

        real(rk),   dimension(:), allocatable   ::  &
            norm_1, norm_2, norm_3, unorm_1, unorm_2, unorm_3, r



        !
        ! Interpolate solution to quadrature nodes
        !
        density         = worker%get_primary_field_element('Density',           'value')
        mom1            = worker%get_primary_field_element('Momentum-1',        'value')
        mom2            = worker%get_primary_field_element('Momentum-2',        'value')
        mom3            = worker%get_primary_field_element('Momentum-3',        'value')
        density_nutilde = worker%get_primary_field_element('Density * NuTilde', 'value')



        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2 = mom2 / worker%coordinate('1','element')
        end if



        !
        ! Compute velocities
        !
        invdensity = ONE/density
        u = mom1*invdensity
        v = mom2*invdensity
        w = mom3*invdensity

        
        !
        ! Compute transport velocity
        !
        r = worker%coordinate('1','element')
        v = v - omega*r

        !
        ! Get normal vector
        !
        norm_1  = worker%normal(1)
        norm_2  = worker%normal(2)
        norm_3  = worker%normal(3)

        unorm_1 = worker%unit_normal(1)
        unorm_2 = worker%unit_normal(2)
        unorm_3 = worker%unit_normal(3)


        !
        ! Compute average flux and field difference.
        ! 
        flux_1 = u*density_nutilde
        flux_2 = v*density_nutilde
        flux_3 = w*density_nutilde


        !
        ! Integrate flux
        !
        call worker%integrate_volume('Density * NuTilde',flux_1,flux_2,flux_3)


    end subroutine compute
    !************************************************************************************************









end module spalart_allmaras_volume_advection
