module rans_bc_diffusion
#include <messenger.h>
    use mod_kinds,              only: ik, rk
    use mod_constants,          only: HALF, ONE, TWO
    use mod_spalart_allmaras,   only: SA_sigma, SA_c_n1
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use mod_rans_efficient
    use DNAD_D
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !------------------------------------------------------------------------------------
    type, public, extends(operator_t)   :: rans_bc_diffusion_t


    contains

        procedure   :: init
        procedure   :: compute

    end type rans_bc_diffusion_t
    !*************************************************************************************






contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rans_bc_diffusion_t),   intent(inout) :: self
        
        ! Set operator name
        call self%set_name('RANS BC Diffusion')

        ! Set operator type
        call self%set_operator_type('BC Diffusive Flux')

        ! Set operator equations
        call self%add_primary_field('Density'   )
        call self%add_primary_field('Momentum-1')
        call self%add_primary_field('Momentum-2')
        call self%add_primary_field('Momentum-3')
        call self%add_primary_field('Energy'    )


    end subroutine init
    !********************************************************************************





    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!  @param[in]      face    face_info_t containing indices on location and misc information
    !!  @param[in]      flux    function_into_t containing info on the function being computed
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(rans_bc_diffusion_t),  intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::                    &
            density, grad1_density, grad2_density, grad3_density,   &
            mom1,    grad1_mom1,    grad2_mom1,    grad3_mom1,      &
            mom2,    grad1_mom2,    grad2_mom2,    grad3_mom2,      &
            mom3,    grad1_mom3,    grad2_mom3,    grad3_mom3,      &
            energy,  grad1_energy,  grad2_energy,  grad3_energy,    &
            u, v, w, invdensity, H, p, T,           &
            grad1_T, grad2_T, grad3_T,              &
            mu,     mu_l,   mu_t,                   &
            lambda, lambda_l,   lambda_t,           &
            k,       k_l,     k_t,                  &
            tau_11, tau_22, tau_33,                 &
            tau_12, tau_13, tau_23,                 &
            flux_1, flux_2, flux_3

        type(AD_D), allocatable, dimension(:)   ::  &
            density_nutilde, grad1_density_nutilde, grad2_density_nutilde, grad3_density_nutilde,   &
            nutilde,         grad1_nutilde,         grad2_nutilde,         grad3_nutilde,           &
            dnutilde_ddensity, dnutilde_ddensitynutilde,                                            &
            f_n1, nu_l, chi, diffusion

        real(rk),   allocatable, dimension(:)   ::  r


        ! Interpolate boundary condition state to face quadrature nodes
        density       = worker%get_field('Density', 'value', 'boundary')
        grad1_density = worker%get_field('Density', 'grad1', 'boundary')
        grad2_density = worker%get_field('Density', 'grad2', 'boundary')
        grad3_density = worker%get_field('Density', 'grad3', 'boundary')

        mom1       = worker%get_field('Momentum-1', 'value', 'boundary')
        grad1_mom1 = worker%get_field('Momentum-1', 'grad1', 'boundary')
        grad2_mom1 = worker%get_field('Momentum-1', 'grad2', 'boundary')
        grad3_mom1 = worker%get_field('Momentum-1', 'grad3', 'boundary')

        mom2       = worker%get_field('Momentum-2', 'value', 'boundary')
        grad1_mom2 = worker%get_field('Momentum-2', 'grad1', 'boundary')
        grad2_mom2 = worker%get_field('Momentum-2', 'grad2', 'boundary')
        grad3_mom2 = worker%get_field('Momentum-2', 'grad3', 'boundary')

        mom3       = worker%get_field('Momentum-3', 'value', 'boundary')
        grad1_mom3 = worker%get_field('Momentum-3', 'grad1', 'boundary')
        grad2_mom3 = worker%get_field('Momentum-3', 'grad2', 'boundary')
        grad3_mom3 = worker%get_field('Momentum-3', 'grad3', 'boundary')

        energy       = worker%get_field('Energy', 'value', 'boundary')
        grad1_energy = worker%get_field('Energy', 'grad1', 'boundary')
        grad2_energy = worker%get_field('Energy', 'grad2', 'boundary')
        grad3_energy = worker%get_field('Energy', 'grad3', 'boundary')

        density_nutilde       = worker%get_field('Density * NuTilde', 'value', 'boundary')
        grad1_density_nutilde = worker%get_field('Density * NuTilde', 'grad1', 'boundary')
        grad2_density_nutilde = worker%get_field('Density * NuTilde', 'grad2', 'boundary')
        grad3_density_nutilde = worker%get_field('Density * NuTilde', 'grad3', 'boundary')

        ! Account for cylindrical. Get tangential momentum from angular momentum.
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','boundary')
            mom2 = mom2 / r
            grad1_mom2 = (grad1_mom2/r) - mom2/r
            grad2_mom2 = (grad2_mom2/r)
            grad3_mom2 = (grad3_mom2/r)
        end if


        !
        ! Compute velocities
        !
        invdensity = ONE/density
        u          = mom1*invdensity
        v          = mom2*invdensity
        w          = mom3*invdensity


        !
        ! Compute pressure
        !
        call compute_pressure_temperature(density,mom1,mom2,mom3,energy,p,T)



        !
        ! Compute viscosity
        !
        call compute_viscosity(density,T,density_nutilde,mu_l, mu_t, lambda_l, lambda_t)
        mu     = mu_l     + mu_t
        lambda = lambda_l + lambda_t




        !
        ! Compute thermal conductivity
        !
        call compute_thermal_conductivity(mu_l,mu_t,k_l,k_t)
        k = k_l + k_t




        !
        ! Compute temperature gradient
        !
        call compute_temperature_gradient(density,u,v,w,p,                              &
                                          grad1_density,grad2_density,grad3_density,    &
                                          grad1_mom1,   grad2_mom1,   grad3_mom1,       &
                                          grad1_mom2,   grad2_mom2,   grad3_mom2,       &
                                          grad1_mom3,   grad2_mom3,   grad3_mom3,       &
                                          grad1_energy, grad2_energy, grad3_energy,     &
                                          grad1_T, grad2_T, grad3_T)


        !
        ! Compute shear stress components
        !
        call compute_shear_stress(worker%coordinate_system(), worker%coordinate('1','boundary'),    &
                                  density,grad1_density,grad2_density,grad3_density,    &
                                  mom1,   grad1_mom1,   grad2_mom1,   grad3_mom1,       &
                                  mom2,   grad1_mom2,   grad2_mom2,   grad3_mom2,       &
                                  mom3,   grad1_mom3,   grad2_mom3,   grad3_mom3,       &
                                  energy, grad1_energy, grad2_energy, grad3_energy,     &
                                  mu, lambda,                                           &
                                  tau_11,tau_22,tau_33,tau_12,tau_13,tau_23)
                                  


        !=================================================
        ! Mass flux
        !=================================================

        !=================================================
        ! momentum-1 flux
        !=================================================
        flux_1 = -tau_11
        flux_2 = -tau_12
        flux_3 = -tau_13
        call worker%integrate_boundary_condition('Momentum-1','Diffusion',flux_1,flux_2,flux_3)

        !=================================================
        ! momentum-2 flux
        !=================================================
        flux_1 = -tau_12
        flux_2 = -tau_22
        flux_3 = -tau_23
        ! Convert to tangential to angular momentum flux
        if (worker%coordinate_system() == 'Cylindrical') then
            flux_1 = flux_1 * r
            flux_2 = flux_2 * r
            flux_3 = flux_3 * r
        end if
        call worker%integrate_boundary_condition('Momentum-2','Diffusion',flux_1,flux_2,flux_3)

        !=================================================
        ! momentum-3 flux
        !=================================================
        flux_1 = -tau_13
        flux_2 = -tau_23
        flux_3 = -tau_33
        call worker%integrate_boundary_condition('Momentum-3','Diffusion',flux_1,flux_2,flux_3)

        !=================================================
        ! Energy flux
        !=================================================
        flux_1 = -k*grad1_T  -  (u*tau_11 + v*tau_12 + w*tau_13)
        flux_2 = -k*grad2_T  -  (u*tau_12 + v*tau_22 + w*tau_23)
        flux_3 = -k*grad3_T  -  (u*tau_13 + v*tau_23 + w*tau_33)
        call worker%integrate_boundary_condition('Energy','Diffusion',flux_1,flux_2,flux_3)



        !=================================================
        ! Spalart-Allmaras: Density * NuTilde
        !=================================================
        ! Compute nutilde, chi
        nu_l    = mu_l*invdensity
        nutilde = density_nutilde*invdensity
        chi     = nutilde/nu_l

        ! Initialize derivatives first
        f_n1 = density
        where (nutilde >= ZERO)
            f_n1 = ONE
        else where
            f_n1 = (SA_c_n1 + chi*chi*chi)/(SA_c_n1 - chi*chi*chi)
        end where


        ! Compute nutilde jacobians for gradient calculation using chain-rule.
        dnutilde_ddensity        = -density_nutilde*invdensity*invdensity
        dnutilde_ddensitynutilde =  invdensity

        grad1_nutilde = dnutilde_ddensity * grad1_density  +  dnutilde_ddensitynutilde * grad1_density_nutilde
        grad2_nutilde = dnutilde_ddensity * grad2_density  +  dnutilde_ddensitynutilde * grad2_density_nutilde
        grad3_nutilde = dnutilde_ddensity * grad3_density  +  dnutilde_ddensitynutilde * grad3_density_nutilde


        diffusion = -(ONE/SA_sigma)*(mu_l + f_n1*density_nutilde)
        flux_1 = diffusion*grad1_nutilde
        flux_2 = diffusion*grad2_nutilde
        flux_3 = diffusion*grad3_nutilde
        call worker%integrate_boundary_condition('Density * NuTilde','Diffusion',flux_1,flux_2,flux_3)



    end subroutine compute
    !**********************************************************************************************









end module rans_bc_diffusion
