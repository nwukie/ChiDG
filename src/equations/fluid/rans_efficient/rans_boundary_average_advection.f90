module rans_boundary_average_advection
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, HALF
    use mod_spalart_allmaras,   only: SA_sigma, SA_c_n1
    use mod_rans_efficient

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none

    private


    !> Implementation of the Euler boundary average flux
    !!
    !!  - At a boundary interface, the solution states Q- and Q+ exists on opposite 
    !!    sides of the boundary. The average flux is computed as Favg = 1/2(F(Q-) + F(Q+))
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/17/2019
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: rans_boundary_average_advection_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rans_boundary_average_advection_t
    !********************************************************************************



contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   1/17/2019
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rans_boundary_average_advection_t),   intent(inout) :: self
        
        ! Set operator name
        call self%set_name('RANS Boundary Average Advection')

        ! Set operator type
        call self%set_operator_type('Boundary Advective Flux')

        ! Set operator equations
        call self%add_primary_field('Density'   )
        call self%add_primary_field('Momentum-1')
        call self%add_primary_field('Momentum-2')
        call self%add_primary_field('Momentum-3')
        call self%add_primary_field('Energy'    )

    end subroutine init
    !********************************************************************************



    !>  Boundary Flux routine for Euler
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/17/2019
    !!
    !!-------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(rans_boundary_average_advection_t),   intent(inout)   :: self
        type(chidg_worker_t),                               intent(inout)   :: worker
        class(properties_t),                                intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::                            &
            density_m, mom1_m, mom2_m, mom3_m, energy_m,                    &
            density_p, mom1_p, mom2_p, mom3_p, energy_p,                    &
            u_m, v_m, w_m, invdensity_m, p_m, temperature_m, H_m,           &
            u_p, v_p, w_p, invdensity_p, p_p, temperature_p, H_p,           &
            flux_1_m, flux_2_m, flux_3_m,                                   &
            flux_1_p, flux_2_p, flux_3_p, r

        type(AD_D), allocatable, dimension(:)   ::  &
            density_nutilde_m, density_nutilde_p



        ! Interpolate solution to quadrature nodes
        density_m = worker%get_field('Density', 'value', 'face interior')
        density_p = worker%get_field('Density', 'value', 'face exterior')

        mom1_m    = worker%get_field('Momentum-1', 'value', 'face interior')
        mom1_p    = worker%get_field('Momentum-1', 'value', 'face exterior')

        mom2_m    = worker%get_field('Momentum-2', 'value', 'face interior')
        mom2_p    = worker%get_field('Momentum-2', 'value', 'face exterior')

        mom3_m    = worker%get_field('Momentum-3', 'value', 'face interior')
        mom3_p    = worker%get_field('Momentum-3', 'value', 'face exterior')

        energy_m  = worker%get_field('Energy', 'value', 'face interior')
        energy_p  = worker%get_field('Energy', 'value', 'face exterior')


        select case (trim(turbulence_model))
            case('Spalart Allmaras')
                density_nutilde_m = worker%get_field('Density * NuTilde', 'value', 'face interior')
                density_nutilde_p = worker%get_field('Density * NuTilde', 'value', 'face exterior')
            case('none')
                density_nutilde_m = ZERO*density_m
                density_nutilde_p = ZERO*density_p
            case default
                call chidg_signal_one(FATAL,"RANS Efficient: invalid turbulence_model.",trim(turbulence_model))
        end select



        ! Account for cylindrical. Get tangential momentum from angular momentum.
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','face interior')
            mom2_m = mom2_m / r
            mom2_p = mom2_p / r
        end if



        !
        ! Compute velocities
        !
        invdensity_m = ONE/density_m
        invdensity_p = ONE/density_p

        u_m = mom1_m*invdensity_m
        v_m = mom2_m*invdensity_m
        w_m = mom3_m*invdensity_m

        u_p = mom1_p*invdensity_p
        v_p = mom2_p*invdensity_p
        w_p = mom3_p*invdensity_p



        !
        ! Compute pressure
        !
        call compute_pressure_temperature(density_m,mom1_m,mom2_m,mom3_m,energy_m,p_m,temperature_m)
        call compute_pressure_temperature(density_p,mom1_p,mom2_p,mom3_p,energy_p,p_p,temperature_p)



        !
        ! Compute boundary condition energy and enthalpy
        !
        H_m = (energy_m + p_m)*invdensity_m
        H_p = (energy_p + p_p)*invdensity_p




        !----------------------------------
        !            mass flux
        !----------------------------------
        flux_1_m = (density_m * u_m)
        flux_2_m = (density_m * v_m)
        flux_3_m = (density_m * w_m)

        flux_1_p = (density_p * u_p)
        flux_2_p = (density_p * v_p)
        flux_3_p = (density_p * w_p)

        call worker%integrate_boundary_average('Density','Advection',       &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)


        !----------------------------------
        !         momentum-1 flux
        !----------------------------------
        flux_1_m = (mom1_m * u_m) + p_m
        flux_2_m = (mom1_m * v_m)
        flux_3_m = (mom1_m * w_m)

        flux_1_p = (mom1_p * u_p) + p_p
        flux_2_p = (mom1_p * v_p)
        flux_3_p = (mom1_p * w_p)

        call worker%integrate_boundary_average('Momentum-1','Advection',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)

        !----------------------------------
        !         momentum-2 flux
        !----------------------------------
        flux_1_m = (mom2_m * u_m)
        flux_2_m = (mom2_m * v_m) + p_m
        flux_3_m = (mom2_m * w_m)
                            
        flux_1_p = (mom2_p * u_p)
        flux_2_p = (mom2_p * v_p) + p_p
        flux_3_p = (mom2_p * w_p)

        ! Convert to tangential to angular momentum flux
        if (worker%coordinate_system() == 'Cylindrical') then
            flux_1_m = flux_1_m * r
            flux_2_m = flux_2_m * r
            flux_3_m = flux_3_m * r

            flux_1_p = flux_1_p * r
            flux_2_p = flux_2_p * r
            flux_3_p = flux_3_p * r
        end if

        call worker%integrate_boundary_average('Momentum-2','Advection',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)


        !----------------------------------
        !         momentum-3 flux
        !----------------------------------
        flux_1_m = (mom3_m * u_m)
        flux_2_m = (mom3_m * v_m)
        flux_3_m = (mom3_m * w_m) + p_m
                    
        flux_1_p = (mom3_p * u_p)
        flux_2_p = (mom3_p * v_p)
        flux_3_p = (mom3_p * w_p) + p_p

        call worker%integrate_boundary_average('Momentum-3','Advection',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)

        !----------------------------------
        !           energy flux
        !----------------------------------
        flux_1_m = (density_m * H_m * u_m)
        flux_2_m = (density_m * H_m * v_m)
        flux_3_m = (density_m * H_m * w_m)

        flux_1_p = (density_p * H_p * u_p)
        flux_2_p = (density_p * H_p * v_p)
        flux_3_p = (density_p * H_p * w_p)

        call worker%integrate_boundary_average('Energy','Advection',        &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)



        !-----------------------------------------
        !            TURBULENCE FLUX
        !-----------------------------------------
        select case (trim(turbulence_model))
            case('Spalart Allmaras')
                flux_1_m = density_nutilde_m * u_m
                flux_2_m = density_nutilde_m * v_m
                flux_3_m = density_nutilde_m * w_m

                flux_1_p = density_nutilde_p * u_p
                flux_2_p = density_nutilde_p * v_p
                flux_3_p = density_nutilde_p * w_p
                call worker%integrate_boundary_average('Density * NuTilde','Advection', &
                                                        flux_1_m, flux_2_m, flux_3_m,   &
                                                        flux_1_p, flux_2_p, flux_3_p)
            case('none')

            case default
                call chidg_signal_one(FATAL,"RANS Efficient: invalid turbulence_model.",trim(turbulence_model))
        end select



    end subroutine compute
    !*********************************************************************************************************












end module rans_boundary_average_advection
