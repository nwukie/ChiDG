module rstm_ssglrrw_volume_advection
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
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: rstm_ssglrrw_volume_advection_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_volume_advection_operator_t
    !***********************************************************************************************

contains


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rstm_ssglrrw_volume_advection_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('RSTMSSGLRRW Volume Advection Operator')

        ! Set operator type
        call self%set_operator_type('Volume Advective Operator')

        ! Set operator equations
        call self%add_primary_field('Density * Omega')
        call self%add_primary_field('Density * Reynolds-11')
        call self%add_primary_field('Density * Reynolds-22')
        call self%add_primary_field('Density * Reynolds-33')
        call self%add_primary_field('Density * Reynolds-12')
        call self%add_primary_field('Density * Reynolds-13')
        call self%add_primary_field('Density * Reynolds-23')



    end subroutine init
    !***********************************************************************************************


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(rstm_ssglrrw_volume_advection_operator_t), intent(inout)   :: self
        type(chidg_worker_t),                                intent(inout)   :: worker
        class(properties_t),                                 intent(inout)   :: prop

        type(AD_D), dimension(:), allocatable   ::          &
            density, mom1, mom2, mom3,  density_omega,    &
            density_reynolds_11, density_reynolds_22, density_reynolds_33, &
            density_reynolds_12, density_reynolds_13, density_reynolds_23, &
            invdensity, u, v, w, flux_1, flux_2, flux_3, r


        !
        ! Interpolate solution to quadrature nodes
        !
        density         = worker%get_field('Density',           'value', 'element')
        mom1            = worker%get_field('Momentum-1',        'value', 'element')
        mom2            = worker%get_field('Momentum-2',        'value', 'element')
        mom3            = worker%get_field('Momentum-3',        'value', 'element')

        density_omega           = worker%get_field('Density * Omega',             'value', 'element')
        density_reynolds_11     = worker%get_field('Density * Reynolds-11',       'value', 'element')
        density_reynolds_22     = worker%get_field('Density * Reynolds-22',       'value', 'element')
        density_reynolds_33     = worker%get_field('Density * Reynolds-33',       'value', 'element')
        density_reynolds_12     = worker%get_field('Density * Reynolds-12',       'value', 'element')
        density_reynolds_13     = worker%get_field('Density * Reynolds-13',       'value', 'element')
        density_reynolds_23     = worker%get_field('Density * Reynolds-23',       'value', 'element')



        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','element') 
            mom2 = mom2 / r
        end if

        
        !
        ! Get fluid advection velocity
        ! 
        invdensity = ONE/density
        u = mom1*invdensity
        v = mom2*invdensity
        w = mom3*invdensity


        !
        ! Compute average flux and field difference.
        ! 
        flux_1 = density_omega * u
        flux_2 = density_omega * v
        flux_3 = density_omega * w


        !
        ! Integrate flux
        !
        call worker%integrate_volume_flux('Density * Omega','Advection',flux_1,flux_2,flux_3)

        !
        ! Compute average flux and field difference.
        ! 
        flux_1 = density_reynolds_11 * u
        flux_2 = density_reynolds_11 * v
        flux_3 = density_reynolds_11 * w


        !
        ! Integrate flux
        !
        call worker%integrate_volume_flux('Density * Reynolds-11','Advection',flux_1,flux_2,flux_3)

        !
        ! Compute average flux and field difference.
        ! 
        flux_1 = density_reynolds_22 * u
        flux_2 = density_reynolds_22 * v
        flux_3 = density_reynolds_22 * w


        !
        ! Integrate flux
        !
        call worker%integrate_volume_flux('Density * Reynolds-22','Advection',flux_1,flux_2,flux_3)

        !
        ! Compute average flux and field difference.
        ! 
        flux_1 = density_reynolds_33 * u
        flux_2 = density_reynolds_33 * v
        flux_3 = density_reynolds_33 * w


        !
        ! Integrate flux
        !
        call worker%integrate_volume_flux('Density * Reynolds-33','Advection',flux_1,flux_2,flux_3)

        !
        ! Compute average flux and field difference.
        ! 
        flux_1 = density_reynolds_12 * u
        flux_2 = density_reynolds_12 * v
        flux_3 = density_reynolds_12 * w


        !
        ! Integrate flux
        !
        call worker%integrate_volume_flux('Density * Reynolds-12','Advection',flux_1,flux_2,flux_3)

        !
        ! Compute average flux and field difference.
        ! 
        flux_1 = density_reynolds_13 * u
        flux_2 = density_reynolds_13 * v
        flux_3 = density_reynolds_13 * w


        !
        ! Integrate flux
        !
        call worker%integrate_volume_flux('Density * Reynolds-13','Advection',flux_1,flux_2,flux_3)

        !
        ! Compute average flux and field difference.
        ! 
        flux_1 = density_reynolds_23 * u
        flux_2 = density_reynolds_23 * v
        flux_3 = density_reynolds_23 * w


        !
        ! Integrate flux
        !
        call worker%integrate_volume_flux('Density * Reynolds-23','Advection',flux_1,flux_2,flux_3)


    end subroutine compute
    !************************************************************************************************









end module rstm_ssglrrw_volume_advection
