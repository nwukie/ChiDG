module rstm_ssglrrw_advection_boundary_average
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
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: rstm_ssglrrw_advection_boundary_average_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_advection_boundary_average_operator_t
    !*****************************************************************************************

contains


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rstm_ssglrrw_advection_boundary_average_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('RSTMSSGLRRW Advection Boundary Average Operator')

        ! Set operator type
        call self%set_operator_type('Boundary Advective Operator')

        ! Set operator equations
        call self%add_primary_field('Density * Omega')
        call self%add_primary_field('Density * Reynolds-11')
        call self%add_primary_field('Density * Reynolds-22')
        call self%add_primary_field('Density * Reynolds-33')
        call self%add_primary_field('Density * Reynolds-12')
        call self%add_primary_field('Density * Reynolds-13')
        call self%add_primary_field('Density * Reynolds-23')


    end subroutine init
    !*****************************************************************************************


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(rstm_ssglrrw_advection_boundary_average_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),                                           intent(inout)   :: worker
        class(properties_t),                                            intent(inout)   :: prop


        type(AD_D), dimension(:), allocatable   ::  &
            density_m, density_p,                   &
            mom1_m, mom2_m, mom3_m,                 &
            mom1_p, mom2_p, mom3_p,                 &
            density_omega_m, density_omega_p,       &
            reynolds_11_m,  reynolds_11_p,          &
            reynolds_22_m,  reynolds_22_p,          &
            reynolds_33_m,  reynolds_33_p,          &
            reynolds_12_m,  reynolds_12_p,          &
            reynolds_13_m,  reynolds_13_p,          &
            reynolds_23_m,  reynolds_23_p,          &
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

        density_omega_m = worker%get_field('Density * Omega', 'value', 'face interior')
        density_omega_p = worker%get_field('Density * Omega', 'value', 'face exterior')

        reynolds_11_m = worker%get_field('Density * Reynolds-11', 'value', 'face interior')
        reynolds_11_p = worker%get_field('Density * Reynolds-11', 'value', 'face exterior')

        reynolds_22_m = worker%get_field('Density * Reynolds-22', 'value', 'face interior')
        reynolds_22_p = worker%get_field('Density * Reynolds-22', 'value', 'face exterior')

        reynolds_33_m = worker%get_field('Density * Reynolds-33', 'value', 'face interior')
        reynolds_33_p = worker%get_field('Density * Reynolds-33', 'value', 'face exterior')
        
        reynolds_12_m = worker%get_field('Density * Reynolds-12', 'value', 'face interior')
        reynolds_12_p = worker%get_field('Density * Reynolds-12', 'value', 'face exterior')

        reynolds_13_m = worker%get_field('Density * Reynolds-13', 'value', 'face interior')
        reynolds_13_p = worker%get_field('Density * Reynolds-13', 'value', 'face exterior')

        reynolds_23_m = worker%get_field('Density * Reynolds-23', 'value', 'face interior')
        reynolds_23_p = worker%get_field('Density * Reynolds-23', 'value', 'face exterior')

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
        flux_1_m = density_omega_m * u_m
        flux_2_m = density_omega_m * v_m
        flux_3_m = density_omega_m * w_m

        flux_1_p = density_omega_p * u_p
        flux_2_p = density_omega_p * v_p
        flux_3_p = density_omega_p * w_p


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_average('Density * Omega','Advection', &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)

        !
        ! Compute average flux
        ! 
        flux_1_m = reynolds_11_m * u_m
        flux_2_m = reynolds_11_m * v_m
        flux_3_m = reynolds_11_m * w_m

        flux_1_p = reynolds_11_p * u_p
        flux_2_p = reynolds_11_p * v_p
        flux_3_p = reynolds_11_p * w_p


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_average('Density * Reynolds-11','Advection', &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)

        !
        ! Compute average flux
        ! 
        flux_1_m = reynolds_22_m * u_m
        flux_2_m = reynolds_22_m * v_m
        flux_3_m = reynolds_22_m * w_m

        flux_1_p = reynolds_22_p * u_p
        flux_2_p = reynolds_22_p * v_p
        flux_3_p = reynolds_22_p * w_p


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_average('Density * Reynolds-22','Advection', &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)

        !
        ! Compute average flux
        ! 
        flux_1_m = reynolds_33_m * u_m
        flux_2_m = reynolds_33_m * v_m
        flux_3_m = reynolds_33_m * w_m

        flux_1_p = reynolds_33_p * u_p
        flux_2_p = reynolds_33_p * v_p
        flux_3_p = reynolds_33_p * w_p


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_average('Density * Reynolds-33','Advection', &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)


        !
        ! Compute average flux
        ! 
        flux_1_m = reynolds_12_m * u_m
        flux_2_m = reynolds_12_m * v_m
        flux_3_m = reynolds_12_m * w_m

        flux_1_p = reynolds_12_p * u_p
        flux_2_p = reynolds_12_p * v_p
        flux_3_p = reynolds_12_p * w_p


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_average('Density * Reynolds-12','Advection', &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)


        !
        ! Compute average flux
        ! 
        flux_1_m = reynolds_13_m * u_m
        flux_2_m = reynolds_13_m * v_m
        flux_3_m = reynolds_13_m * w_m

        flux_1_p = reynolds_13_p * u_p
        flux_2_p = reynolds_13_p * v_p
        flux_3_p = reynolds_13_p * w_p


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_average('Density * Reynolds-13','Advection', &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)

        !
        ! Compute average flux
        ! 
        flux_1_m = reynolds_23_m * u_m
        flux_2_m = reynolds_23_m * v_m
        flux_3_m = reynolds_23_m * w_m

        flux_1_p = reynolds_23_p * u_p
        flux_2_p = reynolds_23_p * v_p
        flux_3_p = reynolds_23_p * w_p


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_average('Density * Reynolds-23','Advection', &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)


       



    end subroutine compute
    !******************************************************************************************









end module rstm_ssglrrw_advection_boundary_average
