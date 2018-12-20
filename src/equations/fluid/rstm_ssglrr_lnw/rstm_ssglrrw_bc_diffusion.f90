module rstm_ssglrrw_bc_diffusion
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D

    use mod_rstm_ssglrrw
    implicit none

    private



    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: rstm_ssglrrw_bc_diffusion_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_bc_diffusion_operator_t
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
        class(rstm_ssglrrw_bc_diffusion_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("RSTMSSGLRRW BC Diffusion Operator")

        !
        ! Set operator type
        !
        call self%set_operator_type("BC Diffusive Flux")

        !
        ! Set operator equations
        !
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
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(rstm_ssglrrw_bc_diffusion_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),                               intent(inout)   :: worker
        class(properties_t),                                intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::            &
            density, mu_l, mu_t, sigma_w,                   &
            diffusion_111, diffusion_112, diffusion_113,            &
            diffusion_221, diffusion_222, diffusion_223,            &
            diffusion_331, diffusion_332, diffusion_333,            &
            diffusion_121, diffusion_122, diffusion_123,            &
            diffusion_131, diffusion_132, diffusion_133,            &
            diffusion_231, diffusion_232, diffusion_233,            &
            grad1_omega, grad2_omega, grad3_omega,          &
            flux_1, flux_2, flux_3, diffusion, integrand

        ! Reynolds stress modelled diffusion flux
        diffusion_111 = worker%get_field('Diffusion-111', 'value', 'boundary')
        diffusion_112 = worker%get_field('Diffusion-112', 'value', 'boundary')
        diffusion_113 = worker%get_field('Diffusion-113', 'value', 'boundary')

        diffusion_221 = worker%get_field('Diffusion-221', 'value', 'boundary')
        diffusion_222 = worker%get_field('Diffusion-222', 'value', 'boundary')
        diffusion_223 = worker%get_field('Diffusion-223', 'value', 'boundary')

        diffusion_331 = worker%get_field('Diffusion-331', 'value', 'boundary')
        diffusion_332 = worker%get_field('Diffusion-332', 'value', 'boundary')
        diffusion_333 = worker%get_field('Diffusion-333', 'value', 'boundary')

        diffusion_121 = worker%get_field('Diffusion-121', 'value', 'boundary')
        diffusion_122 = worker%get_field('Diffusion-122', 'value', 'boundary')
        diffusion_123 = worker%get_field('Diffusion-123', 'value', 'boundary')

        diffusion_131 = worker%get_field('Diffusion-131', 'value', 'boundary')
        diffusion_132 = worker%get_field('Diffusion-132', 'value', 'boundary')
        diffusion_133 = worker%get_field('Diffusion-133', 'value', 'boundary')

        diffusion_231 = worker%get_field('Diffusion-231', 'value', 'boundary')
        diffusion_232 = worker%get_field('Diffusion-232', 'value', 'boundary')
        diffusion_233 = worker%get_field('Diffusion-233', 'value', 'boundary')



        !
        ! Interpolate solution to quadrature nodes
        !
        density             = worker%get_field('Density',          'value','boundary')

        grad1_omega         = worker%get_field('Omega - Gradient 1',          'value','boundary')
        grad2_omega         = worker%get_field('Omega - Gradient 2',          'value','boundary')
        grad3_omega         = worker%get_field('Omega - Gradient 3',          'value','boundary')

        !
        ! Get model fields:
        !   Viscosity
        !
        mu_l = worker%get_field('Laminar Viscosity', 'value', 'boundary')

        mu_t = worker%get_field('Equivalent Eddy Viscosity', 'value', 'boundary')

        sigma_w     = worker%get_field('RSTMSSGLRRW Sigma-w',       'value', 'boundary')
        !-------------------------------------
        !           TURBULENCE FLUX
        !-------------------------------------
        diffusion = -(mu_l + 0.5_rk*mu_t)


        flux_1 = diffusion*grad1_omega
        flux_2 = diffusion*grad2_omega
        flux_3 = diffusion*grad3_omega

        call worker%integrate_boundary_condition('Density * Omega','Diffusion',flux_1,flux_2,flux_3)

        !================================
        !       TURBULENCE FLUX - R_11
        !================================

        flux_1 = -diffusion_111
        flux_2 = -diffusion_112
        flux_3 = -diffusion_113


        call worker%integrate_boundary_condition('Density * Reynolds-11','Diffusion',flux_1,flux_2,flux_3)

        !================================
        !       TURBULENCE FLUX - R_22
        !================================

        flux_1 = -diffusion_221
        flux_2 = -diffusion_222
        flux_3 = -diffusion_223


        call worker%integrate_boundary_condition('Density * Reynolds-22','Diffusion',flux_1,flux_2,flux_3)

        !================================
        !       TURBULENCE FLUX - R_33
        !================================

        flux_1 = -diffusion_331
        flux_2 = -diffusion_332
        flux_3 = -diffusion_333


        call worker%integrate_boundary_condition('Density * Reynolds-33','Diffusion',flux_1,flux_2,flux_3)


        !================================
        !       TURBULENCE FLUX - R_12
        !================================

        flux_1 = -diffusion_121
        flux_2 = -diffusion_122
        flux_3 = -diffusion_123


        call worker%integrate_boundary_condition('Density * Reynolds-12','Diffusion',flux_1,flux_2,flux_3)
        !================================
        !       TURBULENCE FLUX - R_13
        !================================

        flux_1 = -diffusion_131
        flux_2 = -diffusion_132
        flux_3 = -diffusion_133


        call worker%integrate_boundary_condition('Density * Reynolds-13','Diffusion',flux_1,flux_2,flux_3)

        !================================
        !       TURBULENCE FLUX - R_23
        !================================

        flux_1 = -diffusion_231
        flux_2 = -diffusion_232
        flux_3 = -diffusion_233


        call worker%integrate_boundary_condition('Density * Reynolds-23','Diffusion',flux_1,flux_2,flux_3)



    end subroutine compute
    !******************************************************************************************












end module rstm_ssglrrw_bc_diffusion
