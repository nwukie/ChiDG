module rstm_ssglrrw_boundary_diffusion
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF
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
    type, extends(operator_t), public :: rstm_ssglrrw_boundary_diffusion_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_boundary_diffusion_operator_t
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
        class(rstm_ssglrrw_boundary_diffusion_operator_t),   intent(inout) :: self
        
        ! Set operator name
        call self%set_name('RSTMSSGLRRW Boundary Diffusion Operator')

        ! Set operator type
        call self%set_operator_type('Boundary Diffusive Flux')

        ! Set operator equations
        call self%add_primary_field('Density * Omega')
        call self%add_primary_field('Density * Reynolds-11')
        call self%add_primary_field('Density * Reynolds-22')
        call self%add_primary_field('Density * Reynolds-33')
        call self%add_primary_field('Density * Reynolds-12')
        call self%add_primary_field('Density * Reynolds-13')
        call self%add_primary_field('Density * Reynolds-23')

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
        class(rstm_ssglrrw_boundary_diffusion_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),                                   intent(inout)   :: worker
        class(properties_t),                                    intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::            &
            density_m, density_p,                           &
            mu_t_m, mu_l_m,                                 &
            mu_t_p, mu_l_p,                                 &
            sigma_w_m, sigma_w_p,                           &
            grad1_omega_m, grad2_omega_m, grad3_omega_m,    &
            grad1_omega_p, grad2_omega_p, grad3_omega_p,    &
            diffusion_m, diffusion_p,                       &
            diffusion_1_m, diffusion_2_m, diffusion_3_m,    &
            diffusion_1_p, diffusion_2_p, diffusion_3_p,    &
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

        sigma_w_m       = worker%get_field('RSTMSSGLRRW Sigma-w',       'value', 'face interior')
        sigma_w_p       = worker%get_field('RSTMSSGLRRW Sigma-w',       'value', 'face exterior')

        mu_t_m          = worker%get_field('Equivalent Eddy Viscosity', 'value', 'face interior')
        mu_t_p          = worker%get_field('Equivalent Eddy Viscosity', 'value', 'face exterior')

        !-----------------------------------------
        !            TURBULENCE FLUX - Omega
        !-----------------------------------------
        diffusion_m     = (mu_l_m + 0.5_rk*mu_t_m)
        diffusion_p     = (mu_l_p + 0.5_rk*mu_t_p)
        
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
        !            TURBULENCE FLUX - R_11
        !-----------------------------------------
        diffusion_1_m = worker%get_field('Diffusion-111', 'value', 'face interior')
        diffusion_1_p = worker%get_field('Diffusion-111', 'value', 'face exterior')
        
        diffusion_2_m = worker%get_field('Diffusion-112', 'value', 'face interior')
        diffusion_2_p = worker%get_field('Diffusion-112', 'value', 'face exterior')

        diffusion_3_m = worker%get_field('Diffusion-113', 'value', 'face interior')
        diffusion_3_p = worker%get_field('Diffusion-113', 'value', 'face exterior')

        flux_1_m = -diffusion_1_m
        flux_1_p = -diffusion_1_p
        flux_2_m = -diffusion_2_m
        flux_2_p = -diffusion_2_p
        flux_3_m = -diffusion_3_m
        flux_3_p = -diffusion_3_p


        call worker%integrate_boundary_average('Density * Reynolds-11','Diffusion', &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)

        !-----------------------------------------
        !            TURBULENCE FLUX - R_22
        !-----------------------------------------
        diffusion_1_m = worker%get_field('Diffusion-221', 'value', 'face interior')
        diffusion_1_p = worker%get_field('Diffusion-221', 'value', 'face exterior')
        
        diffusion_2_m = worker%get_field('Diffusion-222', 'value', 'face interior')
        diffusion_2_p = worker%get_field('Diffusion-222', 'value', 'face exterior')

        diffusion_3_m = worker%get_field('Diffusion-223', 'value', 'face interior')
        diffusion_3_p = worker%get_field('Diffusion-223', 'value', 'face exterior')

        flux_1_m = -diffusion_1_m
        flux_1_p = -diffusion_1_p
        flux_2_m = -diffusion_2_m
        flux_2_p = -diffusion_2_p
        flux_3_m = -diffusion_3_m
        flux_3_p = -diffusion_3_p


        call worker%integrate_boundary_average('Density * Reynolds-22','Diffusion', &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)

        !-----------------------------------------
        !            TURBULENCE FLUX - R_33
        !-----------------------------------------
        diffusion_1_m = worker%get_field('Diffusion-331', 'value', 'face interior')
        diffusion_1_p = worker%get_field('Diffusion-331', 'value', 'face exterior')
        
        diffusion_2_m = worker%get_field('Diffusion-332', 'value', 'face interior')
        diffusion_2_p = worker%get_field('Diffusion-332', 'value', 'face exterior')

        diffusion_3_m = worker%get_field('Diffusion-333', 'value', 'face interior')
        diffusion_3_p = worker%get_field('Diffusion-333', 'value', 'face exterior')

        flux_1_m = -diffusion_1_m
        flux_1_p = -diffusion_1_p
        flux_2_m = -diffusion_2_m
        flux_2_p = -diffusion_2_p
        flux_3_m = -diffusion_3_m
        flux_3_p = -diffusion_3_p


        call worker%integrate_boundary_average('Density * Reynolds-33','Diffusion', &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)

        !-----------------------------------------
        !            TURBULENCE FLUX - R_12
        !-----------------------------------------
        diffusion_1_m = worker%get_field('Diffusion-121', 'value', 'face interior')
        diffusion_1_p = worker%get_field('Diffusion-121', 'value', 'face exterior')
        
        diffusion_2_m = worker%get_field('Diffusion-122', 'value', 'face interior')
        diffusion_2_p = worker%get_field('Diffusion-122', 'value', 'face exterior')

        diffusion_3_m = worker%get_field('Diffusion-123', 'value', 'face interior')
        diffusion_3_p = worker%get_field('Diffusion-123', 'value', 'face exterior')

        flux_1_m = -diffusion_1_m
        flux_1_p = -diffusion_1_p
        flux_2_m = -diffusion_2_m
        flux_2_p = -diffusion_2_p
        flux_3_m = -diffusion_3_m
        flux_3_p = -diffusion_3_p


        call worker%integrate_boundary_average('Density * Reynolds-12','Diffusion', &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)

        !-----------------------------------------
        !            TURBULENCE FLUX - R_13
        !-----------------------------------------
        diffusion_1_m = worker%get_field('Diffusion-131', 'value', 'face interior')
        diffusion_1_p = worker%get_field('Diffusion-131', 'value', 'face exterior')
        
        diffusion_2_m = worker%get_field('Diffusion-132', 'value', 'face interior')
        diffusion_2_p = worker%get_field('Diffusion-132', 'value', 'face exterior')

        diffusion_3_m = worker%get_field('Diffusion-133', 'value', 'face interior')
        diffusion_3_p = worker%get_field('Diffusion-133', 'value', 'face exterior')

        flux_1_m = -diffusion_1_m
        flux_1_p = -diffusion_1_p
        flux_2_m = -diffusion_2_m
        flux_2_p = -diffusion_2_p
        flux_3_m = -diffusion_3_m
        flux_3_p = -diffusion_3_p


        call worker%integrate_boundary_average('Density * Reynolds-13','Diffusion', &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)

        !-----------------------------------------
        !            TURBULENCE FLUX - R_23
        !-----------------------------------------
        diffusion_1_m = worker%get_field('Diffusion-231', 'value', 'face interior')
        diffusion_1_p = worker%get_field('Diffusion-231', 'value', 'face exterior')
        
        diffusion_2_m = worker%get_field('Diffusion-232', 'value', 'face interior')
        diffusion_2_p = worker%get_field('Diffusion-232', 'value', 'face exterior')

        diffusion_3_m = worker%get_field('Diffusion-233', 'value', 'face interior')
        diffusion_3_p = worker%get_field('Diffusion-233', 'value', 'face exterior')

        flux_1_m = -diffusion_1_m
        flux_1_p = -diffusion_1_p
        flux_2_m = -diffusion_2_m
        flux_2_p = -diffusion_2_p
        flux_3_m = -diffusion_3_m
        flux_3_p = -diffusion_3_p


        call worker%integrate_boundary_average('Density * Reynolds-23','Diffusion', &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)


    end subroutine compute
    !********************************************************************************************












end module rstm_ssglrrw_boundary_diffusion
