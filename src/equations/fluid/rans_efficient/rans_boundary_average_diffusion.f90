module rans_boundary_average_diffusion
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
    type, extends(operator_t), public :: rans_boundary_average_diffusion_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rans_boundary_average_diffusion_t
    !********************************************************************************



contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   1/17/2019
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rans_boundary_average_diffusion_t),   intent(inout) :: self
        
        ! Set operator name
        call self%set_name('RANS Boundary Average Diffusion')

        ! Set operator type
        call self%set_operator_type('Boundary Diffusive Flux')

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
        class(rans_boundary_average_diffusion_t),   intent(inout)   :: self
        type(chidg_worker_t),                               intent(inout)   :: worker
        class(properties_t),                                intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::                            &
            density_m, mom1_m, mom2_m, mom3_m, energy_m,                    &
            density_p, mom1_p, mom2_p, mom3_p, energy_p,                    &
            grad1_density_m, grad2_density_m, grad3_density_m,              &
            grad1_density_p, grad2_density_p, grad3_density_p,              &
            grad1_mom1_m,    grad2_mom1_m,    grad3_mom1_m,                 &
            grad1_mom1_p,    grad2_mom1_p,    grad3_mom1_p,                 &
            grad1_mom2_m,    grad2_mom2_m,    grad3_mom2_m,                 &
            grad1_mom2_p,    grad2_mom2_p,    grad3_mom2_p,                 &
            grad1_mom3_m,    grad2_mom3_m,    grad3_mom3_m,                 &
            grad1_mom3_p,    grad2_mom3_p,    grad3_mom3_p,                 &
            grad1_energy_m,  grad2_energy_m,  grad3_energy_m,               &
            grad1_energy_p,  grad2_energy_p,  grad3_energy_p,               &
            u_m, v_m, w_m, invdensity_m, p_m, temperature_m, H_m,                     &
            u_p, v_p, w_p, invdensity_p, p_p, temperature_p, H_p,                     &
            mu_m, mu_l_m, mu_t_m,                                           &
            mu_p, mu_l_p, mu_t_p,                                           &
            lambda_m, lambda_l_m, lambda_t_m,                               &
            lambda_p, lambda_l_p, lambda_t_p,                               &
            k_m, k_l_m, k_t_m,                                              &
            k_p, k_l_p, k_t_p,                                              &
            grad1_T_m, grad2_T_m, grad3_T_m,                                &
            grad1_T_p, grad2_T_p, grad3_T_p,                                &
            tau_11_m, tau_22_m, tau_33_m, tau_12_m, tau_13_m, tau_23_m,     &
            tau_11_p, tau_22_p, tau_33_p, tau_12_p, tau_13_p, tau_23_p,     &
            flux_1_m, flux_2_m, flux_3_m,                                   &
            flux_1_p, flux_2_p, flux_3_p, r

        type(AD_D), allocatable, dimension(:)   ::  &
            density_nutilde_m, grad1_density_nutilde_m,  grad2_density_nutilde_m,  grad3_density_nutilde_m, &
            density_nutilde_p, grad1_density_nutilde_p,  grad2_density_nutilde_p,  grad3_density_nutilde_p, &
            nutilde_m,         grad1_nutilde_m,         grad2_nutilde_m,         grad3_nutilde_m,           &
            nutilde_p,         grad1_nutilde_p,         grad2_nutilde_p,         grad3_nutilde_p,           &
            dnutilde_ddensity_m, dnutilde_ddensitynutilde_m,                                                &
            dnutilde_ddensity_p, dnutilde_ddensitynutilde_p,                                                &
            f_n1_m, nu_l_m, chi_m, diffusion_m,                                                             &
            f_n1_p, nu_l_p, chi_p, diffusion_p



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



        grad1_density_m    = worker%get_field('Density'   , 'grad1', 'face interior')
        grad2_density_m    = worker%get_field('Density'   , 'grad2', 'face interior')
        grad3_density_m    = worker%get_field('Density'   , 'grad3', 'face interior')

        grad1_mom1_m       = worker%get_field('Momentum-1', 'grad1', 'face interior')
        grad2_mom1_m       = worker%get_field('Momentum-1', 'grad2', 'face interior')
        grad3_mom1_m       = worker%get_field('Momentum-1', 'grad3', 'face interior')

        grad1_mom2_m       = worker%get_field('Momentum-2', 'grad1', 'face interior')
        grad2_mom2_m       = worker%get_field('Momentum-2', 'grad2', 'face interior')
        grad3_mom2_m       = worker%get_field('Momentum-2', 'grad3', 'face interior')

        grad1_mom3_m       = worker%get_field('Momentum-3', 'grad1', 'face interior')
        grad2_mom3_m       = worker%get_field('Momentum-3', 'grad2', 'face interior')
        grad3_mom3_m       = worker%get_field('Momentum-3', 'grad3', 'face interior')

        grad1_energy_m     = worker%get_field('Energy    ', 'grad1', 'face interior')
        grad2_energy_m     = worker%get_field('Energy    ', 'grad2', 'face interior')
        grad3_energy_m     = worker%get_field('Energy    ', 'grad3', 'face interior')





        grad1_density_p    = worker%get_field('Density'   , 'grad1', 'face exterior')
        grad2_density_p    = worker%get_field('Density'   , 'grad2', 'face exterior')
        grad3_density_p    = worker%get_field('Density'   , 'grad3', 'face exterior')

        grad1_mom1_p       = worker%get_field('Momentum-1', 'grad1', 'face exterior')
        grad2_mom1_p       = worker%get_field('Momentum-1', 'grad2', 'face exterior')
        grad3_mom1_p       = worker%get_field('Momentum-1', 'grad3', 'face exterior')

        grad1_mom2_p       = worker%get_field('Momentum-2', 'grad1', 'face exterior')
        grad2_mom2_p       = worker%get_field('Momentum-2', 'grad2', 'face exterior')
        grad3_mom2_p       = worker%get_field('Momentum-2', 'grad3', 'face exterior')

        grad1_mom3_p       = worker%get_field('Momentum-3', 'grad1', 'face exterior')
        grad2_mom3_p       = worker%get_field('Momentum-3', 'grad2', 'face exterior')
        grad3_mom3_p       = worker%get_field('Momentum-3', 'grad3', 'face exterior')

        grad1_energy_p     = worker%get_field('Energy    ', 'grad1', 'face exterior')
        grad2_energy_p     = worker%get_field('Energy    ', 'grad2', 'face exterior')
        grad3_energy_p     = worker%get_field('Energy    ', 'grad3', 'face exterior')



        select case (trim(turbulence_model))
            case('Spalart Allmaras')
                density_nutilde_m       = worker%get_field('Density * NuTilde', 'value', 'face interior')
                grad1_density_nutilde_m = worker%get_field('Density * NuTilde', 'grad1', 'face interior')
                grad2_density_nutilde_m = worker%get_field('Density * NuTilde', 'grad2', 'face interior')
                grad3_density_nutilde_m = worker%get_field('Density * NuTilde', 'grad3', 'face interior')

                density_nutilde_p       = worker%get_field('Density * NuTilde', 'value', 'face exterior')
                grad1_density_nutilde_p = worker%get_field('Density * NuTilde', 'grad1', 'face exterior')
                grad2_density_nutilde_p = worker%get_field('Density * NuTilde', 'grad2', 'face exterior')
                grad3_density_nutilde_p = worker%get_field('Density * NuTilde', 'grad3', 'face exterior')
            case('none')
                density_nutilde_m       = ZERO*density_m 
                grad1_density_nutilde_m = ZERO*density_m 
                grad2_density_nutilde_m = ZERO*density_m 
                grad3_density_nutilde_m = ZERO*density_m 

                density_nutilde_p       = ZERO*density_p 
                grad1_density_nutilde_p = ZERO*density_p 
                grad2_density_nutilde_p = ZERO*density_p 
                grad3_density_nutilde_p = ZERO*density_p 
            case default
                call chidg_signal_one(FATAL,"RANS Efficient: invalid turbulence_model.",trim(turbulence_model))
        end select






        ! Account for cylindrical. Get tangential momentum from angular momentum.
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','face interior')
            mom2_m = mom2_m / r
            mom2_p = mom2_p / r
            grad1_mom2_m = (grad1_mom2_m/r) - mom2_m/r
            grad2_mom2_m = (grad2_mom2_m/r)
            grad3_mom2_m = (grad3_mom2_m/r)
            grad1_mom2_p = (grad1_mom2_p/r) - mom2_p/r
            grad2_mom2_p = (grad2_mom2_p/r)
            grad3_mom2_p = (grad3_mom2_p/r)
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
        ! Compute viscosity
        !
        call compute_viscosity(density_m,temperature_m,density_nutilde_m,mu_l_m, mu_t_m, lambda_l_m, lambda_t_m)
        call compute_viscosity(density_p,temperature_p,density_nutilde_p,mu_l_p, mu_t_p, lambda_l_p, lambda_t_p)
        mu_m     = mu_l_m + mu_t_m
        mu_p     = mu_l_p + mu_t_p
        lambda_m = lambda_l_m + lambda_t_m
        lambda_p = lambda_l_p + lambda_t_p


        !
        ! Compute thermal conductivity
        !
        call compute_thermal_conductivity(mu_l_m,mu_t_m,k_l_m,k_t_m)
        call compute_thermal_conductivity(mu_l_p,mu_t_p,k_l_p,k_t_p)
        k_m = k_l_m + k_t_m
        k_p = k_l_p + k_t_p



        !
        ! Compute temperature gradient
        !
        call compute_temperature_gradient(density_m,u_m,v_m,w_m,p_m,                        &
                                          grad1_density_m,grad2_density_m,grad3_density_m,  &
                                          grad1_mom1_m,   grad2_mom1_m,   grad3_mom1_m,     &
                                          grad1_mom2_m,   grad2_mom2_m,   grad3_mom2_m,     &
                                          grad1_mom3_m,   grad2_mom3_m,   grad3_mom3_m,     &
                                          grad1_energy_m, grad2_energy_m, grad3_energy_m,   &
                                          grad1_T_m, grad2_T_m, grad3_T_m)

        call compute_temperature_gradient(density_p,u_p,v_p,w_p,p_p,                        &
                                          grad1_density_p,grad2_density_p,grad3_density_p,  &
                                          grad1_mom1_p,   grad2_mom1_p,   grad3_mom1_p,     &
                                          grad1_mom2_p,   grad2_mom2_p,   grad3_mom2_p,     &
                                          grad1_mom3_p,   grad2_mom3_p,   grad3_mom3_p,     &
                                          grad1_energy_p, grad2_energy_p, grad3_energy_p,   &
                                          grad1_T_p, grad2_T_p, grad3_T_p)



        !
        ! Compute shear stress components
        !
        call compute_shear_stress(worker%coordinate_system(), worker%coordinate('1','face interior'),   &
                                  density_m,grad1_density_m,grad2_density_m,grad3_density_m,            &
                                  mom1_m,   grad1_mom1_m,   grad2_mom1_m,   grad3_mom1_m,               &
                                  mom2_m,   grad1_mom2_m,   grad2_mom2_m,   grad3_mom2_m,               &
                                  mom3_m,   grad1_mom3_m,   grad2_mom3_m,   grad3_mom3_m,               &
                                  energy_m, grad1_energy_m, grad2_energy_m, grad3_energy_m,             &
                                  mu_m, lambda_m,                                                       &
                                  tau_11_m,tau_22_m,tau_33_m,tau_12_m,tau_13_m,tau_23_m)
                                  
        call compute_shear_stress(worker%coordinate_system(), worker%coordinate('1','face interior'),   &
                                  density_p,grad1_density_p,grad2_density_p,grad3_density_p,            &
                                  mom1_p,   grad1_mom1_p,   grad2_mom1_p,   grad3_mom1_p,               &
                                  mom2_p,   grad1_mom2_p,   grad2_mom2_p,   grad3_mom2_p,               &
                                  mom3_p,   grad1_mom3_p,   grad2_mom3_p,   grad3_mom3_p,               &
                                  energy_p, grad1_energy_p, grad2_energy_p, grad3_energy_p,             &
                                  mu_p, lambda_p,                                                       &
                                  tau_11_p,tau_22_p,tau_33_p,tau_12_p,tau_13_p,tau_23_p)





        !----------------------------------
        !            mass flux
        !----------------------------------




        !----------------------------------
        !         momentum-1 flux
        !----------------------------------
        flux_1_m = -tau_11_m
        flux_2_m = -tau_12_m
        flux_3_m = -tau_13_m

        flux_1_p = -tau_11_p
        flux_2_p = -tau_12_p
        flux_3_p = -tau_13_p

        call worker%integrate_boundary_average('Momentum-1','Diffusion',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)


        !----------------------------------
        !         momentum-2 flux
        !----------------------------------
        flux_1_m = -tau_12_m
        flux_2_m = -tau_22_m
        flux_3_m = -tau_23_m

        flux_1_p = -tau_12_p
        flux_2_p = -tau_22_p
        flux_3_p = -tau_23_p

        ! Convert to tangential to angular momentum flux
        if (worker%coordinate_system() == 'Cylindrical') then
            flux_1_m = flux_1_m * r
            flux_2_m = flux_2_m * r
            flux_3_m = flux_3_m * r

            flux_1_p = flux_1_p * r
            flux_2_p = flux_2_p * r
            flux_3_p = flux_3_p * r
        end if

        call worker%integrate_boundary_average('Momentum-2','Diffusion',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)

        !----------------------------------
        !         momentum-3 flux
        !----------------------------------
        flux_1_m = -tau_13_m
        flux_2_m = -tau_23_m
        flux_3_m = -tau_33_m

        flux_1_p = -tau_13_p
        flux_2_p = -tau_23_p
        flux_3_p = -tau_33_p

        call worker%integrate_boundary_average('Momentum-3','Diffusion',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)

        !----------------------------------
        !           energy flux
        !----------------------------------
        flux_1_m = -k_m*grad1_T_m  -  (u_m*tau_11_m + v_m*tau_12_m + w_m*tau_13_m)
        flux_2_m = -k_m*grad2_T_m  -  (u_m*tau_12_m + v_m*tau_22_m + w_m*tau_23_m)
        flux_3_m = -k_m*grad3_T_m  -  (u_m*tau_13_m + v_m*tau_23_m + w_m*tau_33_m)

        flux_1_p = -k_p*grad1_T_p  -  (u_p*tau_11_p + v_p*tau_12_p + w_p*tau_13_p)
        flux_2_p = -k_p*grad2_T_p  -  (u_p*tau_12_p + v_p*tau_22_p + w_p*tau_23_p)
        flux_3_p = -k_p*grad3_T_p  -  (u_p*tau_13_p + v_p*tau_23_p + w_p*tau_33_p)

        call worker%integrate_boundary_average('Energy','Diffusion',        &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)




        !-----------------------------------------
        !            TURBULENCE FLUX
        !-----------------------------------------
        select case (trim(turbulence_model))
            case('Spalart Allmaras')
                nu_l_m = mu_l_m*invdensity_m
                nu_l_p = mu_l_p*invdensity_p

                nutilde_m = density_nutilde_m*invdensity_m
                nutilde_p = density_nutilde_p*invdensity_p

                chi_m = nutilde_m/nu_l_m
                chi_p = nutilde_p/nu_l_p


                ! Initialize derivatives first
                f_n1_m = density_m
                where (nutilde_m >= ZERO)
                    f_n1_m = ONE
                else where
                    f_n1_m = (SA_c_n1 + chi_m*chi_m*chi_m)/(SA_c_n1 - chi_m*chi_m*chi_m)
                end where

                f_n1_p = density_p
                where (nutilde_p >= ZERO)
                    f_n1_p = ONE
                else where
                    f_n1_p = (SA_c_n1 + chi_p*chi_p*chi_p)/(SA_c_n1 - chi_p*chi_p*chi_p)
                end where


                ! Compute nutilde jacobians for gradient calculation using chain-rule.
                dnutilde_ddensity_m        = -density_nutilde_m*invdensity_m*invdensity_m
                dnutilde_ddensitynutilde_m =  invdensity_m

                dnutilde_ddensity_p        = -density_nutilde_p*invdensity_p*invdensity_p
                dnutilde_ddensitynutilde_p =  invdensity_p

                grad1_nutilde_m = dnutilde_ddensity_m * grad1_density_m  +  dnutilde_ddensitynutilde_m * grad1_density_nutilde_m
                grad2_nutilde_m = dnutilde_ddensity_m * grad2_density_m  +  dnutilde_ddensitynutilde_m * grad2_density_nutilde_m
                grad3_nutilde_m = dnutilde_ddensity_m * grad3_density_m  +  dnutilde_ddensitynutilde_m * grad3_density_nutilde_m

                grad1_nutilde_p = dnutilde_ddensity_p * grad1_density_p  +  dnutilde_ddensitynutilde_p * grad1_density_nutilde_p
                grad2_nutilde_p = dnutilde_ddensity_p * grad2_density_p  +  dnutilde_ddensitynutilde_p * grad2_density_nutilde_p
                grad3_nutilde_p = dnutilde_ddensity_p * grad3_density_p  +  dnutilde_ddensitynutilde_p * grad3_density_nutilde_p


                ! DIFFUSION
                diffusion_m = -(ONE/SA_sigma)*(mu_l_m + f_n1_m*density_nutilde_m)
                diffusion_p = -(ONE/SA_sigma)*(mu_l_p + f_n1_p*density_nutilde_p)

                flux_1_m = diffusion_m*grad1_nutilde_m
                flux_2_m = diffusion_m*grad2_nutilde_m
                flux_3_m = diffusion_m*grad3_nutilde_m

                flux_1_p = diffusion_p*grad1_nutilde_p
                flux_2_p = diffusion_p*grad2_nutilde_p
                flux_3_p = diffusion_p*grad3_nutilde_p
                call worker%integrate_boundary_average('Density * NuTilde','Diffusion', &
                                                        flux_1_m, flux_2_m, flux_3_m,   &
                                                        flux_1_p, flux_2_p, flux_3_p)

            case('none')

            case default
                call chidg_signal_one(FATAL,"RANS Efficient: invalid turbulence_model.",trim(turbulence_model))
        end select




    end subroutine compute
    !*********************************************************************************************************












end module rans_boundary_average_diffusion
