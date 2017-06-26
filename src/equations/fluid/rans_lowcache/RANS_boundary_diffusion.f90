module RANS_boundary_diffusion
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, HALF, ZERO, THREE, FOUR
    use mod_fluid,              only: omega
    use mod_spalart_allmaras

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
    !!  @date   1/28/2016
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: RANS_boundary_diffusion_t

        real(rk)    :: gam = 1.4_rk
        real(rk)    :: R   = 287.15_rk
        real(rk)    :: Cp  = 1003.0_rk
        real(rk)    :: Pr  = 0.72_rk

    contains

        procedure   :: init
        procedure   :: compute

    end type RANS_boundary_diffusion_t
    !********************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(RANS_boundary_diffusion_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('RANS Boundary Diffusion')

        !
        ! Set operator type
        !
        call self%set_operator_type('Boundary Diffusive Flux')

        !
        ! Set operator equations
        !
        call self%add_primary_field('Density'          )
        call self%add_primary_field('Momentum-1'       )
        call self%add_primary_field('Momentum-2'       )
        call self%add_primary_field('Momentum-3'       )
        call self%add_primary_field('Energy'           )
        call self%add_primary_field('Density * NuTilde')

    end subroutine init
    !********************************************************************************



    !>  Boundary Flux routine for Euler
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!-------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(RANS_boundary_diffusion_t),   intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::                                        &
            density_m, mom1_m, mom2_m, mom3_m, energy_m, density_nutilde_m,             &
            density_p, mom1_p, mom2_p, mom3_p, energy_p, density_nutilde_p,             &
            grad1_density_m,         grad2_density_m,         grad3_density_m,          &
            grad1_density_p,         grad2_density_p,         grad3_density_p,          &
            grad1_mom1_m,            grad2_mom1_m,            grad3_mom1_m,             &
            grad1_mom2_m,            grad2_mom2_m,            grad3_mom2_m,             &
            grad1_mom3_m,            grad2_mom3_m,            grad3_mom3_m,             &
            grad1_mom1_p,            grad2_mom1_p,            grad3_mom1_p,             &
            grad1_mom2_p,            grad2_mom2_p,            grad3_mom2_p,             &
            grad1_mom3_p,            grad2_mom3_p,            grad3_mom3_p,             &
            grad1_energy_m,          grad2_energy_m,          grad3_energy_m,           &
            grad1_energy_p,          grad2_energy_p,          grad3_energy_p,           &
            grad1_density_nutilde_m, grad2_density_nutilde_m, grad3_density_nutilde_m,  &
            grad1_density_nutilde_p, grad2_density_nutilde_p, grad3_density_nutilde_p,  &
            grad1_nutilde_m,         grad2_nutilde_m,         grad3_nutilde_m,          &
            grad1_nutilde_p,         grad2_nutilde_p,         grad3_nutilde_p,          &
            grad1_T_m, grad2_T_m, grad3_T_m,                                    &
            grad1_T_p, grad2_T_p, grad3_T_p,                                    &
            grad1_u_m, grad2_u_m, grad3_u_m,                                    &
            grad1_v_m, grad2_v_m, grad3_v_m,                                    &
            grad1_w_m, grad2_w_m, grad3_w_m,                                    &
            grad1_u_p, grad2_u_p, grad3_u_p,                                    &
            grad1_v_p, grad2_v_p, grad3_v_p,                                    &
            grad1_w_p, grad2_w_p, grad3_w_p,                                    &
            dnutilde_ddensity_m,        dnutilde_ddensity_p,                    &
            dnutilde_ddensitynutilde_m, dnutilde_ddensitynutilde_p,             &
            dke_ddensity_m, dke_dmom1_m, dke_dmom2_m, dke_dmom3_m,              &
            dke_ddensity_p, dke_dmom1_p, dke_dmom2_p, dke_dmom3_p,              &
            dp_ddensity_m, dp_dmom1_m, dp_dmom2_m, dp_dmom3_m, dp_denergy_m,    &
            dp_ddensity_p, dp_dmom1_p, dp_dmom2_p, dp_dmom3_p, dp_denergy_p,    &
            dT_ddensity_m, dT_dmom1_m, dT_dmom2_m, dT_dmom3_m, dT_denergy_m,    &
            dT_ddensity_p, dT_dmom1_p, dT_dmom2_p, dT_dmom3_p, dT_denergy_p,    &
            du_ddensity_m, dv_ddensity_m, dw_ddensity_m,                        &
            du_ddensity_p, dv_ddensity_p, dw_ddensity_p,                        &
            du_dmom1_m, dv_dmom2_m, dw_dmom3_m,                                 &
            du_dmom1_p, dv_dmom2_p, dw_dmom3_p,                                 &
            vorticity_1_m, vorticity_2_m, vorticity_3_m,                        &
            vorticity_1_p, vorticity_2_p, vorticity_3_p,                        &
            div_velocity_m, div_velocity_p,                                     &
            u_m, v_m, w_m, invdensity_m,                                        &
            u_p, v_p, w_p, invdensity_p,                                        &
            p_m, p_p,                                                           &
            mu_l_m,  mu_l_p,  mu_t_m,  mu_t_p,  mu_m,  mu_p,                    &
            mu2_l_m, mu2_l_p, mu2_t_m, mu2_t_p, mu2_m, mu2_p,                   &
            k_l_m,   k_l_p,   k_t_m,   k_t_p,   k_m,   k_p,                     &
            chi_m, chi_p, f_n1_m, f_n1_p, f_v1_m, f_v1_p,                       &
            nu_l_m, nu_l_p, nutilde_m, nutilde_p,                               &
            tau_11_m, tau_22_m, tau_33_m, tau_12_m, tau_13_m, tau_23_m,         &
            tau_11_p, tau_22_p, tau_33_p, tau_12_p, tau_13_p, tau_23_p,         &
            diffusion_m, diffusion_p,                                           &
            flux_1_m, flux_2_m, flux_3_m,                                       &
            flux_1_p, flux_2_p, flux_3_p,                                       &
            flux_1, flux_2, flux_3, integrand


        real(rk), allocatable, dimension(:) ::      &
            norm_1, norm_2, norm_3, r

        real(rk)                            :: const

        ! Sutherlands Law constants
        real(rk) :: mu0 = 1.7894e-5_rk  ! [kg/(m*s)]
        real(rk) :: T0  = 273.11_rk     ! [K]
        real(rk) :: S   = 110.56_rk     ! [K]



        !
        ! Interpolate solution to quadrature nodes
        !
        density_m         = worker%get_primary_field_face('Density'   , 'value', 'face interior')
        density_p         = worker%get_primary_field_face('Density'   , 'value', 'face exterior')

        mom1_m            = worker%get_primary_field_face('Momentum-1', 'value', 'face interior')
        mom1_p            = worker%get_primary_field_face('Momentum-1', 'value', 'face exterior')

        mom2_m            = worker%get_primary_field_face('Momentum-2', 'value', 'face interior')
        mom2_p            = worker%get_primary_field_face('Momentum-2', 'value', 'face exterior')

        mom3_m            = worker%get_primary_field_face('Momentum-3', 'value', 'face interior')
        mom3_p            = worker%get_primary_field_face('Momentum-3', 'value', 'face exterior')

        energy_m          = worker%get_primary_field_face('Energy'    , 'value', 'face interior')
        energy_p          = worker%get_primary_field_face('Energy'    , 'value', 'face exterior')

        density_nutilde_m = worker%get_primary_field_face('Density * NuTilde', 'value', 'face interior')
        density_nutilde_p = worker%get_primary_field_face('Density * NuTilde', 'value', 'face exterior')
        invdensity_m = ONE/density_m
        invdensity_p = ONE/density_p



        !
        ! Interpolate gradient to quadrature nodes: Interior face
        !
        grad1_density_m         = worker%get_primary_field_face('Density'   ,       'grad1+lift', 'face interior')
        grad2_density_m         = worker%get_primary_field_face('Density'   ,       'grad2+lift', 'face interior')
        grad3_density_m         = worker%get_primary_field_face('Density'   ,       'grad3+lift', 'face interior')

        grad1_mom1_m            = worker%get_primary_field_face('Momentum-1',       'grad1+lift', 'face interior')
        grad2_mom1_m            = worker%get_primary_field_face('Momentum-1',       'grad2+lift', 'face interior')
        grad3_mom1_m            = worker%get_primary_field_face('Momentum-1',       'grad3+lift', 'face interior')

        grad1_mom2_m            = worker%get_primary_field_face('Momentum-2',       'grad1+lift', 'face interior')
        grad2_mom2_m            = worker%get_primary_field_face('Momentum-2',       'grad2+lift', 'face interior')
        grad3_mom2_m            = worker%get_primary_field_face('Momentum-2',       'grad3+lift', 'face interior')

        grad1_mom3_m            = worker%get_primary_field_face('Momentum-3',       'grad1+lift', 'face interior')
        grad2_mom3_m            = worker%get_primary_field_face('Momentum-3',       'grad2+lift', 'face interior')
        grad3_mom3_m            = worker%get_primary_field_face('Momentum-3',       'grad3+lift', 'face interior')

        grad1_energy_m          = worker%get_primary_field_face('Energy'    ,       'grad1+lift', 'face interior')
        grad2_energy_m          = worker%get_primary_field_face('Energy'    ,       'grad2+lift', 'face interior')
        grad3_energy_m          = worker%get_primary_field_face('Energy'    ,       'grad3+lift', 'face interior')

        grad1_density_nutilde_m = worker%get_primary_field_face('Density * NuTilde','grad1+lift', 'face interior')
        grad2_density_nutilde_m = worker%get_primary_field_face('Density * NuTilde','grad2+lift', 'face interior')
        grad3_density_nutilde_m = worker%get_primary_field_face('Density * NuTilde','grad3+lift', 'face interior')


        !
        ! Interpolate gradient to quadrature nodes: Exterior face
        !
        grad1_density_p         = worker%get_primary_field_face('Density'   ,       'grad1+lift', 'face exterior')
        grad2_density_p         = worker%get_primary_field_face('Density'   ,       'grad2+lift', 'face exterior')
        grad3_density_p         = worker%get_primary_field_face('Density'   ,       'grad3+lift', 'face exterior')

        grad1_mom1_p            = worker%get_primary_field_face('Momentum-1',       'grad1+lift', 'face exterior')
        grad2_mom1_p            = worker%get_primary_field_face('Momentum-1',       'grad2+lift', 'face exterior')
        grad3_mom1_p            = worker%get_primary_field_face('Momentum-1',       'grad3+lift', 'face exterior')

        grad1_mom2_p            = worker%get_primary_field_face('Momentum-2',       'grad1+lift', 'face exterior')
        grad2_mom2_p            = worker%get_primary_field_face('Momentum-2',       'grad2+lift', 'face exterior')
        grad3_mom2_p            = worker%get_primary_field_face('Momentum-2',       'grad3+lift', 'face exterior')

        grad1_mom3_p            = worker%get_primary_field_face('Momentum-3',       'grad1+lift', 'face exterior')
        grad2_mom3_p            = worker%get_primary_field_face('Momentum-3',       'grad2+lift', 'face exterior')
        grad3_mom3_p            = worker%get_primary_field_face('Momentum-3',       'grad3+lift', 'face exterior')

        grad1_energy_p          = worker%get_primary_field_face('Energy'    ,       'grad1+lift', 'face exterior')
        grad2_energy_p          = worker%get_primary_field_face('Energy'    ,       'grad2+lift', 'face exterior')
        grad3_energy_p          = worker%get_primary_field_face('Energy'    ,       'grad3+lift', 'face exterior')

        grad1_density_nutilde_p = worker%get_primary_field_face('Density * NuTilde','grad1+lift', 'face exterior')
        grad2_density_nutilde_p = worker%get_primary_field_face('Density * NuTilde','grad2+lift', 'face exterior')
        grad3_density_nutilde_p = worker%get_primary_field_face('Density * NuTilde','grad3+lift', 'face exterior')






        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        r = worker%coordinate('1','boundary')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / r
            mom2_p = mom2_p / r

            grad1_mom2_m = (grad1_mom2_m/r) - mom2_m/r
            grad2_mom2_m = (grad2_mom2_m/r)
            grad3_mom2_m = (grad3_mom2_m/r)

            grad1_mom2_p = (grad1_mom2_p/r) - mom2_p/r
            grad2_mom2_p = (grad2_mom2_p/r)
            grad3_mom2_p = (grad3_mom2_p/r)

        else if (worker%coordinate_system() == 'Cartesian') then

        else
            call chidg_signal(FATAL,"inlet, bad coordinate system")
        end if




        !
        ! Get normal vector
        !
        norm_1 = worker%normal(1)
        norm_2 = worker%normal(2)
        norm_3 = worker%normal(3)




        !
        ! Compute Pressure, Temperature
        !
        p_m = worker%get_model_field_face('Pressure', 'value', 'face interior')
        p_p = worker%get_model_field_face('Pressure', 'value', 'face exterior')


        !
        ! Compute material viscosity, second viscosity, and thermal conductivity
        !
        mu_l_m = worker%get_model_field_face('Laminar Viscosity', 'value', 'face interior')
        mu_l_p = worker%get_model_field_face('Laminar Viscosity', 'value', 'face exterior')
        mu2_l_m = -(TWO/THREE)*mu_l_m
        mu2_l_p = -(TWO/THREE)*mu_l_p
        k_l_m   = self%Cp * mu_l_m / self%Pr
        k_l_p   = self%Cp * mu_l_p / self%Pr



        !
        ! Compute turbulence fields
        !
        nu_l_m    = mu_l_m*invdensity_m
        nu_l_p    = mu_l_p*invdensity_p
        nutilde_m = density_nutilde_m*invdensity_m
        nutilde_p = density_nutilde_p*invdensity_p


        !
        ! Compute chi, f_v1
        !
        chi_m = nutilde_m/nu_l_m
        chi_p = nutilde_p/nu_l_p
        f_v1_m = chi_m*chi_m*chi_m/(chi_m*chi_m*chi_m + SA_c_v1*SA_c_v1*SA_c_v1)
        f_v1_p = chi_p*chi_p*chi_p/(chi_p*chi_p*chi_p + SA_c_v1*SA_c_v1*SA_c_v1)


        !
        ! Compute f_n1_m
        !
        f_n1_m = density_m
        where (nutilde_m >= ZERO)
            f_n1_m = ONE
        else where
            f_n1_m = (SA_c_n1 + chi_m*chi_m*chi_m)/(SA_c_n1 - chi_m*chi_m*chi_m)
        end where

        !
        ! Compute f_n1_p
        !
        f_n1_p = density_p
        where (nutilde_p >= ZERO)
            f_n1_p = ONE
        else where
            f_n1_p = (SA_c_n1 + chi_p*chi_p*chi_p)/(SA_c_n1 - chi_p*chi_p*chi_p)
        end where





        !
        ! Initialize derivatives, compute mu_t
        !
        mu_t_m = density_m
        where (nutilde_m >= 0)
            mu_t_m = density_m * nutilde_m * f_v1_m
        else where
            mu_t_m = ZERO
        end where

        mu_t_p = density_p
        where (nutilde_p >= 0)
            mu_t_p = density_p * nutilde_p * f_v1_p
        else where
            mu_t_p = ZERO
        end where


        !
        ! Compute: 
        !   - Second Coefficient of Turbulent Viscosity, Stokes' Hypothesis.
        !   - Turbulent Thermal Conductivity, Reynolds' analogy.
        !
        mu2_t_m = (-TWO/THREE)*mu_t_m
        mu2_t_p = (-TWO/THREE)*mu_t_p
        k_t_m   = self%Cp*mu_t_m/SA_Pr_t
        k_t_p   = self%Cp*mu_t_p/SA_Pr_t



        !
        ! Compute total coefficients
        !
        mu_m  = mu_l_m  + mu_t_m
        mu_p  = mu_l_p  + mu_t_p
        mu2_m = mu2_l_m + mu2_t_m
        mu2_p = mu2_l_p + mu2_t_p
        k_m   = k_l_m   + k_t_m
        k_p   = k_l_p   + k_t_p



        !
        ! Compute nutilde jacobians for gradient calculation using chain-rule.
        !
        dnutilde_ddensity_m        = -density_nutilde_m*invdensity_m*invdensity_m
        dnutilde_ddensity_p        = -density_nutilde_p*invdensity_p*invdensity_p
        dnutilde_ddensitynutilde_m =  invdensity_m
        dnutilde_ddensitynutilde_p =  invdensity_p


        grad1_nutilde_m = dnutilde_ddensity_m * grad1_density_m  +  dnutilde_ddensitynutilde_m * grad1_density_nutilde_m
        grad2_nutilde_m = dnutilde_ddensity_m * grad2_density_m  +  dnutilde_ddensitynutilde_m * grad2_density_nutilde_m
        grad3_nutilde_m = dnutilde_ddensity_m * grad3_density_m  +  dnutilde_ddensitynutilde_m * grad3_density_nutilde_m

        grad1_nutilde_p = dnutilde_ddensity_p * grad1_density_p  +  dnutilde_ddensitynutilde_p * grad1_density_nutilde_p
        grad2_nutilde_p = dnutilde_ddensity_p * grad2_density_p  +  dnutilde_ddensitynutilde_p * grad2_density_nutilde_p
        grad3_nutilde_p = dnutilde_ddensity_p * grad3_density_p  +  dnutilde_ddensitynutilde_p * grad3_density_nutilde_p



        !############### Computing Temperature Gradient ###################

        !
        ! Compute temperature gradient
        !

        ! Compute velocities
        u_m = mom1_m*invdensity_m
        v_m = mom2_m*invdensity_m
        w_m = mom3_m*invdensity_m

        u_p = mom1_p*invdensity_p
        v_p = mom2_p*invdensity_p
        w_p = mom3_p*invdensity_p

        ! Compute Kinetic Energy Jacobians
        dke_ddensity_m = -HALF*(u_m*u_m + v_m*v_m + w_m*w_m)
        dke_dmom1_m    = u_m
        dke_dmom2_m    = v_m
        dke_dmom3_m    = w_m

        dke_ddensity_p = -HALF*(u_p*u_p + v_p*v_p + w_p*w_p)
        dke_dmom1_p    = u_p
        dke_dmom2_p    = v_p
        dke_dmom3_p    = w_p

        ! Compute Pressure Jacobians
        dp_ddensity_m = -(self%gam-ONE)*dke_ddensity_m
        dp_dmom1_m    = -(self%gam-ONE)*dke_dmom1_m
        dp_dmom2_m    = -(self%gam-ONE)*dke_dmom2_m
        dp_dmom3_m    = -(self%gam-ONE)*dke_dmom3_m
        dp_denergy_m  =  dp_dmom3_m    ! Initialize derivatives
        dp_denergy_m  =  (self%gam-ONE)   ! No negative sign

        dp_ddensity_p = -(self%gam-ONE)*dke_ddensity_p
        dp_dmom1_p    = -(self%gam-ONE)*dke_dmom1_p
        dp_dmom2_p    = -(self%gam-ONE)*dke_dmom2_p
        dp_dmom3_p    = -(self%gam-ONE)*dke_dmom3_p
        dp_denergy_p  =  dp_dmom3_p    ! Initialize derivatives
        dp_denergy_p  =  (self%gam-ONE)   ! No negative sign

        ! Compute Temperature Jacobians
        const = ONE/self%R
        dT_ddensity_m = const*invdensity_m*dp_ddensity_m  -  const*invdensity_m*invdensity_m*p_m
        dT_dmom1_m    = const*invdensity_m*dp_dmom1_m
        dT_dmom2_m    = const*invdensity_m*dp_dmom2_m
        dT_dmom3_m    = const*invdensity_m*dp_dmom3_m
        dT_denergy_m  = const*invdensity_m*dp_denergy_m

        dT_ddensity_p = const*invdensity_p*dp_ddensity_p  -  const*invdensity_p*invdensity_p*p_p
        dT_dmom1_p    = const*invdensity_p*dp_dmom1_p
        dT_dmom2_p    = const*invdensity_p*dp_dmom2_p
        dT_dmom3_p    = const*invdensity_p*dp_dmom3_p
        dT_denergy_p  = const*invdensity_p*dp_denergy_p


        ! Temperature gradient
        grad1_T_m = dT_ddensity_m*grad1_density_m + dT_dmom1_m*grad1_mom1_m + dT_dmom2_m*grad1_mom2_m + dT_dmom3_m*grad1_mom3_m + dT_denergy_m*grad1_energy_m
        grad2_T_m = dT_ddensity_m*grad2_density_m + dT_dmom1_m*grad2_mom1_m + dT_dmom2_m*grad2_mom2_m + dT_dmom3_m*grad2_mom3_m + dT_denergy_m*grad2_energy_m
        grad3_T_m = dT_ddensity_m*grad3_density_m + dT_dmom1_m*grad3_mom1_m + dT_dmom2_m*grad3_mom2_m + dT_dmom3_m*grad3_mom3_m + dT_denergy_m*grad3_energy_m

        grad1_T_p = dT_ddensity_p*grad1_density_p + dT_dmom1_p*grad1_mom1_p + dT_dmom2_p*grad1_mom2_p + dT_dmom3_p*grad1_mom3_p + dT_denergy_p*grad1_energy_p
        grad2_T_p = dT_ddensity_p*grad2_density_p + dT_dmom1_p*grad2_mom1_p + dT_dmom2_p*grad2_mom2_p + dT_dmom3_p*grad2_mom3_p + dT_denergy_p*grad2_energy_p
        grad3_T_p = dT_ddensity_p*grad3_density_p + dT_dmom1_p*grad3_mom1_p + dT_dmom2_p*grad3_mom2_p + dT_dmom3_p*grad3_mom3_p + dT_denergy_p*grad3_energy_p



        !############### Computing Vorticity ######################

        !
        ! compute velocity jacobians
        !
        du_ddensity_m = -invdensity_m*invdensity_m*mom1_m
        dv_ddensity_m = -invdensity_m*invdensity_m*mom2_m
        dw_ddensity_m = -invdensity_m*invdensity_m*mom3_m

        du_ddensity_p = -invdensity_p*invdensity_p*mom1_p
        dv_ddensity_p = -invdensity_p*invdensity_p*mom2_p
        dw_ddensity_p = -invdensity_p*invdensity_p*mom3_p


        du_dmom1_m = invdensity_m
        dv_dmom2_m = invdensity_m
        dw_dmom3_m = invdensity_m

        du_dmom1_p = invdensity_p
        dv_dmom2_p = invdensity_p
        dw_dmom3_p = invdensity_p


        !
        ! compute velocity gradients via chain rule:
        !
        !   u = f(rho,rhou)
        !
        !   grad(u) = dudrho * grad(rho)  +  dudrhou * grad(rhou)
        !
        grad1_u_m = du_ddensity_m*grad1_density_m  +  du_dmom1_m*grad1_mom1_m
        grad2_u_m = du_ddensity_m*grad2_density_m  +  du_dmom1_m*grad2_mom1_m
        grad3_u_m = du_ddensity_m*grad3_density_m  +  du_dmom1_m*grad3_mom1_m

        grad1_v_m = dv_ddensity_m*grad1_density_m  +  dv_dmom2_m*grad1_mom2_m
        grad2_v_m = dv_ddensity_m*grad2_density_m  +  dv_dmom2_m*grad2_mom2_m
        grad3_v_m = dv_ddensity_m*grad3_density_m  +  dv_dmom2_m*grad3_mom2_m

        grad1_w_m = dw_ddensity_m*grad1_density_m  +  dw_dmom3_m*grad1_mom3_m
        grad2_w_m = dw_ddensity_m*grad2_density_m  +  dw_dmom3_m*grad2_mom3_m
        grad3_w_m = dw_ddensity_m*grad3_density_m  +  dw_dmom3_m*grad3_mom3_m


        grad1_u_p = du_ddensity_p*grad1_density_p  +  du_dmom1_p*grad1_mom1_p
        grad2_u_p = du_ddensity_p*grad2_density_p  +  du_dmom1_p*grad2_mom1_p
        grad3_u_p = du_ddensity_p*grad3_density_p  +  du_dmom1_p*grad3_mom1_p

        grad1_v_p = dv_ddensity_p*grad1_density_p  +  dv_dmom2_p*grad1_mom2_p
        grad2_v_p = dv_ddensity_p*grad2_density_p  +  dv_dmom2_p*grad2_mom2_p
        grad3_v_p = dv_ddensity_p*grad3_density_p  +  dv_dmom2_p*grad3_mom2_p

        grad1_w_p = dw_ddensity_p*grad1_density_p  +  dw_dmom3_p*grad1_mom3_p
        grad2_w_p = dw_ddensity_p*grad2_density_p  +  dw_dmom3_p*grad2_mom3_p
        grad3_w_p = dw_ddensity_p*grad3_density_p  +  dw_dmom3_p*grad3_mom3_p




        !----------------------------------------------------------
        !
        !                        Cartesian
        !
        !----------------------------------------------------------
        if (worker%coordinate_system() == 'Cartesian') then



            !
            ! Compute vorticity:
            !
            !   vorticity = Curl(U)
            !
            !   U = [u_x, u_y, u_z] = [u,v,w] 
            !
            !   curl(U) = (dwdx - dvdz)i  +  (dudz - dwdx)j  +  (dvdx - dudy)k
            !
            vorticity_1_m =  (grad2_w_m - grad3_v_m)
            vorticity_2_m =  (grad3_u_m - grad1_w_m)
            vorticity_3_m =  (grad1_v_m - grad2_u_m) 

            vorticity_1_p =  (grad2_w_p - grad3_v_p)
            vorticity_2_p =  (grad3_u_p - grad1_w_p)
            vorticity_3_p =  (grad1_v_p - grad2_u_p) 


        else if (worker%coordinate_system() == 'Cylindrical') then

            !
            ! Compute divergence of velocity vector
            !
            !   U = [u_r, u_theta, u_z] = [u,v,w] 
            !
            !   curl(U) = ((1/r)dwdtheta - dvdz)i  +  (dudz - dwdr)j  +  (1/r)(d(rv)dr - dudtheta)k
            !           = ((1/r)dwdtheta - dvdz)i  +  (dudz - dwdr)j  +  ( dvdr - (1/r)dudtheta  + (v/r) )k
            !           = (grad2_w - grad3_v)i  +  (grad3_u - grad1_w)j  +  ( grad1_v - grad2_u  + (v/r) )k
            !
            !   Note:
            !       grad1_ = d/dr
            !       grad2_ = (1/r)d/dtheta
            !       grad3_ = d/dz
            !
            v_m = mom2_m*invdensity_m
            v_p = mom2_p*invdensity_p

            vorticity_1_m =  (grad2_w_m - grad3_v_m)
            vorticity_2_m =  (grad3_u_m - grad1_w_m)
            vorticity_3_m =  (grad1_v_m - grad2_u_m + (v_m/r)) 

            vorticity_1_p =  (grad2_w_p - grad3_v_p)
            vorticity_2_p =  (grad3_u_p - grad1_w_p)
            vorticity_3_p =  (grad1_v_p - grad2_u_p + (v_p/r)) 



            !
            ! Account for rotation, convert to relative vorticity
            !
            vorticity_3_m = vorticity_3_m - TWO*omega
            vorticity_3_p = vorticity_3_p - TWO*omega

        end if






        !################## Computing Shear Stresses #####################

        !----------------------------------------------------------
        !
        !                        Cartesian
        !
        !----------------------------------------------------------
        if (worker%coordinate_system() == 'Cartesian') then



            !
            ! Compute divergence of velocity vector
            !
            !   U = [u_x, u_y, u_z] = [u,v,w] 
            !
            !   div(V) = dudx + dvdy + dwdz
            !
            div_velocity_m = grad1_u_m + grad2_v_m + grad3_w_m
            div_velocity_p = grad1_u_p + grad2_v_p + grad3_w_p

            !
            ! Compute shear stress components
            !
            tau_11_m = TWO*mu_m*grad1_u_m  +  mu2_m*(div_velocity_m)
            tau_22_m = TWO*mu_m*grad2_v_m  +  mu2_m*(div_velocity_m)
            tau_33_m = TWO*mu_m*grad3_w_m  +  mu2_m*(div_velocity_m)

            tau_12_m = mu_m*(grad2_u_m + grad1_v_m)
            tau_13_m = mu_m*(grad3_u_m + grad1_w_m)
            tau_23_m = mu_m*(grad2_w_m + grad3_v_m)


            tau_11_p = TWO*mu_p*grad1_u_p  +  mu2_p*(div_velocity_p)
            tau_22_p = TWO*mu_p*grad2_v_p  +  mu2_p*(div_velocity_p)
            tau_33_p = TWO*mu_p*grad3_w_p  +  mu2_p*(div_velocity_p)

            tau_12_p = mu_p*(grad2_u_p + grad1_v_p)
            tau_13_p = mu_p*(grad3_u_p + grad1_w_p)
            tau_23_p = mu_p*(grad2_w_p + grad3_v_p)



        else if (worker%coordinate_system() == 'Cylindrical') then

            !
            ! Compute divergence of velocity vector
            !
            !   U = [u_r, u_theta, u_z] = [u,v,w] 
            !
            !   div(U) = (1/r)d(r*u)dr   + (1/r)dv/dtheta + dw/dz
            !          = (du/dr + u/r) + (1/r)dv/dtheta + dw/dz
            !
            !   Note:
            !       grad1_u = du/dr
            !       grad2_v = (1/r)dv/dtheta
            !       grad3_w = dw/dz
            !
            u_m = mom1_m*invdensity_m
            u_p = mom1_p*invdensity_p
            v_m = mom2_m*invdensity_m
            v_p = mom2_p*invdensity_p
            div_velocity_m = (grad1_u_m + u_m/r) + grad2_v_m + grad3_w_m
            div_velocity_p = (grad1_u_p + u_p/r) + grad2_v_p + grad3_w_p

            !
            ! Compute shear stress components
            !
            tau_11_m = TWO*mu_m*(grad1_u_m          )  +  mu2_m*(div_velocity_m)
            tau_22_m = TWO*mu_m*(grad2_v_m + (u_m/r))  +  mu2_m*(div_velocity_m)
            tau_33_m = TWO*mu_m*(grad3_w_m          )  +  mu2_m*(div_velocity_m)

            tau_12_m = mu_m*(grad2_u_m + grad1_v_m  - (v_m/r))
            tau_13_m = mu_m*(grad3_u_m + grad1_w_m           )
            tau_23_m = mu_m*(grad2_w_m + grad3_v_m           )


            tau_11_p = TWO*mu_p*(grad1_u_p          )  +  mu2_p*(div_velocity_p)
            tau_22_p = TWO*mu_p*(grad2_v_p + (u_p/r))  +  mu2_p*(div_velocity_p)
            tau_33_p = TWO*mu_p*(grad3_w_p          )  +  mu2_p*(div_velocity_p)

            tau_12_p = mu_p*(grad2_u_p + grad1_v_p  - (v_p/r))
            tau_13_p = mu_p*(grad3_u_p + grad1_w_p           )
            tau_23_p = mu_p*(grad2_w_p + grad3_v_p           )

        end if























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

        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)


        ! dot with normal vector
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)

        call worker%integrate_boundary('Momentum-1',integrand)


        !----------------------------------
        !         momentum-2 flux
        !----------------------------------
        flux_1_m = -tau_12_m
        flux_2_m = -tau_22_m
        flux_3_m = -tau_23_m

        flux_1_p = -tau_12_p
        flux_2_p = -tau_22_p
        flux_3_p = -tau_23_p

        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)


        ! dot with normal vector
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)

        !
        ! Convert to tangential to angular momentum flux
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            integrand = integrand * r
        else if (worker%coordinate_system() == 'Cartesian') then

        else
            call chidg_signal(FATAL,"inlet, bad coordinate system")
        end if


        call worker%integrate_boundary('Momentum-2',integrand)


        !----------------------------------
        !         momentum-3 flux
        !----------------------------------
        flux_1_m = -tau_13_m
        flux_2_m = -tau_23_m
        flux_3_m = -tau_33_m

        flux_1_p = -tau_13_p
        flux_2_p = -tau_23_p
        flux_3_p = -tau_33_p


        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)


        ! dot with normal vector
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)

        call worker%integrate_boundary('Momentum-3',integrand)


        !----------------------------------
        !           energy flux
        !----------------------------------
        flux_1_m = -k_m*grad1_T_m  -  (u_m*tau_11_m + v_m*tau_12_m + w_m*tau_13_m)
        flux_2_m = -k_m*grad2_T_m  -  (u_m*tau_12_m + v_m*tau_22_m + w_m*tau_23_m)
        flux_3_m = -k_m*grad3_T_m  -  (u_m*tau_13_m + v_m*tau_23_m + w_m*tau_33_m)

        flux_1_p = -k_p*grad1_T_p  -  (u_p*tau_11_p + v_p*tau_12_p + w_p*tau_13_p)
        flux_2_p = -k_p*grad2_T_p  -  (u_p*tau_12_p + v_p*tau_22_p + w_p*tau_23_p)
        flux_3_p = -k_p*grad3_T_p  -  (u_p*tau_13_p + v_p*tau_23_p + w_p*tau_33_p)


        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)


        ! dot with normal vector
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)

        call worker%integrate_boundary('Energy',integrand)


        !-----------------------------------------
        !           turbulence flux
        !-----------------------------------------
        diffusion_m = -(ONE/SA_sigma)*(mu_l_m + f_n1_m*density_nutilde_m)
        diffusion_p = -(ONE/SA_sigma)*(mu_l_p + f_n1_p*density_nutilde_p)

        flux_1_m = diffusion_m*grad1_nutilde_m
        flux_2_m = diffusion_m*grad2_nutilde_m
        flux_3_m = diffusion_m*grad3_nutilde_m

        flux_1_p = diffusion_p*grad1_nutilde_p
        flux_2_p = diffusion_p*grad2_nutilde_p
        flux_3_p = diffusion_p*grad3_nutilde_p


        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)


        ! dot with normal vector
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)

        call worker%integrate_boundary('Density * NuTilde',integrand)







    end subroutine compute
    !*********************************************************************************************************












end module RANS_boundary_diffusion
