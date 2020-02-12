module spalart_allmaras_advection_boundary_average
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF
    use mod_fluid,              only: omega, gam, Rgas
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
    !-----------------------------------------------------------------------------------------
    type, extends(operator_t), public :: spalart_allmaras_advection_boundary_average_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type spalart_allmaras_advection_boundary_average_operator_t
    !*****************************************************************************************

contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init(self)
        class(spalart_allmaras_advection_boundary_average_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('Spalart-Allmaras Advection Boundary Average Operator')

        ! Set operator type
        call self%set_operator_type('Boundary Advective Operator')

        ! Set operator equations
        call self%add_primary_field('Density * NuTilde')

    end subroutine init
    !*****************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(spalart_allmaras_advection_boundary_average_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),                                           intent(inout)   :: worker
        class(properties_t),                                            intent(inout)   :: prop


        type(AD_D), dimension(:), allocatable   ::  &
            density_m, density_p,                   &
            mom1_m, mom2_m, mom3_m,                 &
            mom1_p, mom2_p, mom3_p,                 &
            density_nutilde_m, density_nutilde_p,   &
            invdensity_m, invdensity_p,             &
            u_m, v_m, w_m,                          &
            u_p, v_p, w_p,                          &
            flux_1_m, flux_2_m, flux_3_m,           &
            flux_1_p, flux_2_p, flux_3_p, r


        !
        ! Interpolate solution to quadrature nodes
        !
        density_m         = worker%get_field('Density',           'value', 'face interior')
        density_p         = worker%get_field('Density',           'value', 'face exterior')

        mom1_m            = worker%get_field('Momentum-1',        'value', 'face interior')
        mom1_p            = worker%get_field('Momentum-1',        'value', 'face exterior')

        mom2_m            = worker%get_field('Momentum-2',        'value', 'face interior')
        mom2_p            = worker%get_field('Momentum-2',        'value', 'face exterior')

        mom3_m            = worker%get_field('Momentum-3',        'value', 'face interior')
        mom3_p            = worker%get_field('Momentum-3',        'value', 'face exterior')

        density_nutilde_m = worker%get_field('Density * NuTilde', 'value', 'face interior')
        density_nutilde_p = worker%get_field('Density * NuTilde', 'value', 'face exterior')


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','face interior') 
            mom2_m = mom2_m / r
            mom2_p = mom2_p / r
        end if



        !
        ! Get fluid advection velocity
        !
        invdensity_m = ONE/density_m
        u_m = mom1_m*invdensity_m
        v_m = mom2_m*invdensity_m
        w_m = mom3_m*invdensity_m

        invdensity_p = ONE/density_p
        u_p = mom1_p*invdensity_p
        v_p = mom2_p*invdensity_p
        w_p = mom3_p*invdensity_p


        !
        ! Compute average flux
        ! 
        flux_1_m = density_nutilde_m * u_m
        flux_2_m = density_nutilde_m * v_m
        flux_3_m = density_nutilde_m * w_m

        flux_1_p = density_nutilde_p * u_p
        flux_2_p = density_nutilde_p * v_p
        flux_3_p = density_nutilde_p * w_p


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_average('Density * NuTilde','Advection', &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)


    end subroutine compute
    !******************************************************************************************









end module spalart_allmaras_advection_boundary_average
