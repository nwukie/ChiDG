module sst_boundary_diffusion
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use mod_fluid
    use DNAD_D
    implicit none

    private



    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: sst_boundary_diffusion_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type sst_boundary_diffusion_operator_t
    !********************************************************************************


contains



    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(sst_boundary_diffusion_operator_t),   intent(inout) :: self
        
        ! Set operator name
        call self%set_name('SST Boundary Diffusion Operator')

        ! Set operator type
        call self%set_operator_type('Boundary Diffusive Flux')

        ! Set operator equations
        call self%add_primary_field('Density * Omega')
        call self%add_primary_field('Density * k')
        call self%add_primary_field('Energy')
        
    end subroutine init
    !********************************************************************************



    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(sst_boundary_diffusion_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),                                   intent(inout)   :: worker
        class(properties_t),                                    intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::            &
            density_m, density_p,                           &
            mu_t_m, mu_l_m,                                 &
            mu_t_p, mu_l_p,                                 &
            sigma_w_m, sigma_w_p,                           &
            sigma_k_m, sigma_k_p,                           &
            grad1_omega_m, grad2_omega_m, grad3_omega_m,    &
            grad1_omega_p, grad2_omega_p, grad3_omega_p,    &
            grad1_k_m, grad2_k_m, grad3_k_m,    &
            grad1_k_p, grad2_k_p, grad3_k_p,    &
            grad1_temp_m, grad2_temp_m, grad3_temp_m,    &
            grad1_temp_p, grad2_temp_p, grad3_temp_p,    &
            diffusion_m, diffusion_p,                       &
            diffusion_1_m, diffusion_2_m, diffusion_3_m,    &
            diffusion_1_p, diffusion_2_p, diffusion_3_p,    &
            energy_flux_1_m, energy_flux_2_m, energy_flux_3_m,    &
            energy_flux_1_p, energy_flux_2_p, energy_flux_3_p,    &
            flux_1_m, flux_2_m, flux_3_m,                   &
            flux_1_p, flux_2_p, flux_3_p

        !
        ! Interpolate solution to quadrature nodes
        !
        density_m   = worker%get_field('Density',               'value', 'face interior')
        density_p   = worker%get_field('Density',               'value', 'face exterior')


        mu_l_m      = worker%get_field('Laminar Viscosity',     'value', 'face interior')
        mu_l_p      = worker%get_field('Laminar Viscosity',     'value', 'face exterior')


        ! Omega diffusive terms

        grad1_omega_m   = worker%get_field('Omega - Gradient 1',        'value', 'face interior')
        grad1_omega_p   = worker%get_field('Omega - Gradient 1',        'value', 'face exterior')
        
        grad2_omega_m   = worker%get_field('Omega - Gradient 2',        'value', 'face interior')
        grad2_omega_p   = worker%get_field('Omega - Gradient 2',        'value', 'face exterior')

        grad3_omega_m   = worker%get_field('Omega - Gradient 3',        'value', 'face interior')
        grad3_omega_p   = worker%get_field('Omega - Gradient 3',        'value', 'face exterior')

        grad1_k_m   = worker%get_field('k - Gradient 1',        'value', 'face interior')
        grad1_k_p   = worker%get_field('k - Gradient 1',        'value', 'face exterior')
        
        grad2_k_m   = worker%get_field('k - Gradient 2',        'value', 'face interior')
        grad2_k_p   = worker%get_field('k - Gradient 2',        'value', 'face exterior')

        grad3_k_m   = worker%get_field('k - Gradient 3',        'value', 'face interior')
        grad3_k_p   = worker%get_field('k - Gradient 3',        'value', 'face exterior')


        sigma_w_m       = worker%get_field('SST sigma_w',       'value', 'face interior')
        sigma_w_p       = worker%get_field('SST sigma_w',       'value', 'face exterior')

        sigma_k_m       = worker%get_field('SST sigma_k',       'value', 'face interior')
        sigma_k_p       = worker%get_field('SST sigma_k',       'value', 'face exterior')

        mu_t_m          = worker%get_field('Turbulent Viscosity', 'value', 'face interior')
        mu_t_p          = worker%get_field('Turbulent Viscosity', 'value', 'face exterior')

!        grad1_temp_m   = worker%get_field('Temperature Gradient - 1',        'value', 'face interior')
!        grad1_temp_p   = worker%get_field('Temperature Gradient - 1',        'value', 'face exterior')
!        
!        grad2_temp_m   = worker%get_field('Temperature Gradient - 2',        'value', 'face interior')
!        grad2_temp_p   = worker%get_field('Temperature Gradient - 2',        'value', 'face exterior')
!
!        grad3_temp_m   = worker%get_field('Temperature Gradient - 3',        'value', 'face interior')
!        grad3_temp_p   = worker%get_field('Temperature Gradient - 3',        'value', 'face exterior')
!
!
!        energy_flux_1_m =  cp*mu_t_m*grad1_temp_m/0.9_rk
!        energy_flux_2_m =  cp*mu_t_m*grad2_temp_m/0.9_rk
!        energy_flux_3_m =  cp*mu_t_m*grad3_temp_m/0.9_rk
!
!
!        energy_flux_1_p =  cp*mu_t_p*grad1_temp_p/0.9_rk
!        energy_flux_2_p =  cp*mu_t_p*grad2_temp_p/0.9_rk
!        energy_flux_3_p =  cp*mu_t_p*grad3_temp_p/0.9_rk

        !-----------------------------------------
        !            TURBULENCE FLUX - Omega
        !-----------------------------------------
        diffusion_m     = (mu_l_m + sigma_w_m*mu_t_m)
        diffusion_p     = (mu_l_p + sigma_w_p*mu_t_p)
        
        flux_1_m = -diffusion_m*grad1_omega_m
        flux_1_p = -diffusion_p*grad1_omega_p
        flux_2_m = -diffusion_m*grad2_omega_m
        flux_2_p = -diffusion_p*grad2_omega_p
        flux_3_m = -diffusion_m*grad3_omega_m
        flux_3_p = -diffusion_p*grad3_omega_p

        call worker%integrate_boundary_average('Density * Omega','Diffusion', &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)

        !-----------------------------------------
        !            TURBULENCE FLUX - k
        !-----------------------------------------
        diffusion_m     = (mu_l_m + sigma_k_m*mu_t_m)
        diffusion_p     = (mu_l_p + sigma_k_p*mu_t_p)
        
        flux_1_m = -diffusion_m*grad1_k_m
        flux_1_p = -diffusion_p*grad1_k_p
        flux_2_m = -diffusion_m*grad2_k_m
        flux_2_p = -diffusion_p*grad2_k_p
        flux_3_m = -diffusion_m*grad3_k_m
        flux_3_p = -diffusion_p*grad3_k_p

        call worker%integrate_boundary_average('Density * k','Diffusion', &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)



        !================================
        !       TURBULENCE FLUX - Energy 
        !================================

        !flux_1_m = -energy_flux_1_m
        !flux_2_m = -energy_flux_2_m
        !flux_3_m = -energy_flux_3_m

        !flux_1_p = -energy_flux_1_p
        !flux_2_p = -energy_flux_2_p
        !flux_3_p = -energy_flux_3_p

        flux_1_m = -diffusion_m*grad1_k_m
        flux_1_p = -diffusion_p*grad1_k_p
        flux_2_m = -diffusion_m*grad2_k_m
        flux_2_p = -diffusion_p*grad2_k_p
        flux_3_m = -diffusion_m*grad3_k_m
        flux_3_p = -diffusion_p*grad3_k_p



        call worker%integrate_boundary_average('Energy','Diffusion', &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)


    end subroutine compute
    !********************************************************************************************












end module sst_boundary_diffusion
