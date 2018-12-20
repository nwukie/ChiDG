module sst_bc_diffusion
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF
    use mod_fluid

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
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
    type, extends(operator_t), public :: sst_bc_diffusion_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type sst_bc_diffusion_operator_t
    !*****************************************************************************************










contains



    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(sst_bc_diffusion_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("SST BC Diffusion Operator")

        !
        ! Set operator type
        !
        call self%set_operator_type("BC Diffusive Flux")

        !
        ! Set operator equations
        !
        call self%add_primary_field('Density * Omega')
        call self%add_primary_field('Density * k')
        call self%add_primary_field('Energy')


    end subroutine init
    !*****************************************************************************************



    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(sst_bc_diffusion_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),                               intent(inout)   :: worker
        class(properties_t),                                intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::            &
            density, mu_l, mu_t, sigma_w, sigma_k,                   &
            grad1_omega, grad2_omega, grad3_omega,          &
            grad1_k, grad2_k, grad3_k,          &
            grad1_temp, grad2_temp, grad3_temp,          &
            energy_flux_1, energy_flux_2, energy_flux_3, &
            flux_1, flux_2, flux_3, diffusion, integrand


        !
        ! Interpolate solution to quadrature nodes
        !
        density             = worker%get_field('Density',          'value','boundary')

        grad1_omega         = worker%get_field('Omega - Gradient 1',          'value','boundary')
        grad2_omega         = worker%get_field('Omega - Gradient 2',          'value','boundary')
        grad3_omega         = worker%get_field('Omega - Gradient 3',          'value','boundary')

        grad1_k         = worker%get_field('k - Gradient 1',          'value','boundary')
        grad2_k         = worker%get_field('k - Gradient 2',          'value','boundary')
        grad3_k         = worker%get_field('k - Gradient 3',          'value','boundary')


        !
        ! Get model fields:
        !   Viscosity
        !
        mu_l = worker%get_field('Laminar Viscosity', 'value', 'boundary')

        mu_t = worker%get_field('Turbulent Viscosity', 'value', 'boundary')

        sigma_w     = worker%get_field('SST sigma_w',       'value', 'boundary')
        sigma_k     = worker%get_field('SST sigma_k',       'value', 'boundary')


        !grad1_temp = worker%get_field('Temperature Gradient - 1', 'value')
        !grad2_temp = worker%get_field('Temperature Gradient - 2', 'value')
        !grad3_temp = worker%get_field('Temperature Gradient - 3', 'value')


        !energy_flux_1 = HALF*(diffusion_111+diffusion_221+diffusion_331) + cp*mu_t*grad1_temp/0.9_rk
        !energy_flux_2 = HALF*(diffusion_112+diffusion_222+diffusion_332) + cp*mu_t*grad2_temp/0.9_rk
        !energy_flux_3 = HALF*(diffusion_113+diffusion_223+diffusion_333) + cp*mu_t*grad3_temp/0.9_rk


        !-------------------------------------
        !           TURBULENCE FLUX
        !-------------------------------------
        diffusion = -(mu_l + sigma_w*mu_t)


        flux_1 = diffusion*grad1_omega
        flux_2 = diffusion*grad2_omega
        flux_3 = diffusion*grad3_omega

        call worker%integrate_boundary_condition('Density * Omega','Diffusion',flux_1,flux_2,flux_3)

        !-------------------------------------
        !           TURBULENCE FLUX
        !-------------------------------------
        diffusion = -(mu_l + sigma_k*mu_t)


        flux_1 = diffusion*grad1_k
        flux_2 = diffusion*grad2_k
        flux_3 = diffusion*grad3_k

        call worker%integrate_boundary_condition('Density * k','Diffusion',flux_1,flux_2,flux_3)


        !================================
        !       TURBULENCE FLUX - Energy 
        !================================

        !flux_1 = -energy_flux_1
        !flux_2 = -energy_flux_2
        !flux_3 = -energy_flux_3

        flux_1 = diffusion*grad1_k
        flux_2 = diffusion*grad2_k
        flux_3 = diffusion*grad3_k


        call worker%integrate_boundary_condition('Energy','Diffusion',flux_1,flux_2,flux_3)



    end subroutine compute
    !******************************************************************************************












end module sst_bc_diffusion
