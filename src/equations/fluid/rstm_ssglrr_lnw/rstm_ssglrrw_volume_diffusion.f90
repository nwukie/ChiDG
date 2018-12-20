module rstm_ssglrrw_volume_diffusion
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
    !! @date    01/30/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: rstm_ssglrrw_volume_diffusion_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_volume_diffusion_operator_t
    !********************************************************************************










contains


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/30/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rstm_ssglrrw_volume_diffusion_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('RSTMSSGLRRW Volume Diffusion Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('Volume Diffusive Flux')

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
        !call self%add_primary_field('Energy')

    end subroutine init
    !********************************************************************************


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/30/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(rstm_ssglrrw_volume_diffusion_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),                                   intent(inout)   :: worker
        class(properties_t),                                    intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::                    &
            density, invdensity,      mu_l,                         &
            diffusion_111, diffusion_112, diffusion_113,            &
            diffusion_221, diffusion_222, diffusion_223,            &
            diffusion_331, diffusion_332, diffusion_333,            &
            diffusion_121, diffusion_122, diffusion_123,            &
            diffusion_131, diffusion_132, diffusion_133,            &
            diffusion_231, diffusion_232, diffusion_233,            &
            grad1_omega, grad2_omega, grad3_omega,                  &
            grad1_temp, grad2_temp, grad3_temp,                  &
            energy_flux_1, energy_flux_2, energy_flux_3,            &
            sigma_w, mu_t,                                          &
            flux_1, flux_2, flux_3, diffusion


        !
        ! Get model fields:
        !   Viscosity
        !
        mu_l = worker%get_field('Laminar Viscosity', 'value', 'element')


        sigma_w     = worker%get_field('RSTMSSGLRRW Sigma-w',       'value')

        mu_t        = worker%get_field('Equivalent Eddy Viscosity', 'value')

        ! Reynolds stress modelled diffusion flux
        diffusion_111 = worker%get_field('Diffusion-111', 'value')
        diffusion_112 = worker%get_field('Diffusion-112', 'value')
        diffusion_113 = worker%get_field('Diffusion-113', 'value')

        diffusion_221 = worker%get_field('Diffusion-221', 'value')
        diffusion_222 = worker%get_field('Diffusion-222', 'value')
        diffusion_223 = worker%get_field('Diffusion-223', 'value')

        diffusion_331 = worker%get_field('Diffusion-331', 'value')
        diffusion_332 = worker%get_field('Diffusion-332', 'value')
        diffusion_333 = worker%get_field('Diffusion-333', 'value')

        diffusion_121 = worker%get_field('Diffusion-121', 'value')
        diffusion_122 = worker%get_field('Diffusion-122', 'value')
        diffusion_123 = worker%get_field('Diffusion-123', 'value')

        diffusion_131 = worker%get_field('Diffusion-131', 'value')
        diffusion_132 = worker%get_field('Diffusion-132', 'value')
        diffusion_133 = worker%get_field('Diffusion-133', 'value')

        diffusion_231 = worker%get_field('Diffusion-231', 'value')
        diffusion_232 = worker%get_field('Diffusion-232', 'value')
        diffusion_233 = worker%get_field('Diffusion-233', 'value')

        ! Turbulent heat flux for the energy equation
        ! is modeled as one half of the trace of the diffusion tensor
        grad1_temp = worker%get_field('Temperature Gradient - 1', 'value')
        grad2_temp = worker%get_field('Temperature Gradient - 2', 'value')
        grad3_temp = worker%get_field('Temperature Gradient - 3', 'value')


        energy_flux_1 = HALF*(diffusion_111+diffusion_221+diffusion_331) + cp*mu_t*grad1_temp/0.9_rk
        energy_flux_2 = HALF*(diffusion_112+diffusion_222+diffusion_332) + cp*mu_t*grad2_temp/0.9_rk
        energy_flux_3 = HALF*(diffusion_113+diffusion_223+diffusion_333) + cp*mu_t*grad3_temp/0.9_rk

        ! Omega diffusive terms

        grad1_omega = worker%get_field('Omega - Gradient 1',        'value')
        grad2_omega = worker%get_field('Omega - Gradient 2',        'value')
        grad3_omega = worker%get_field('Omega - Gradient 3',        'value')

                


        !================================
        !       TURBULENCE FLUX - Omega
        !================================
        diffusion = (mu_l + 0.5_rk*mu_t)

        flux_1 = -diffusion*grad1_omega
        flux_2 = -diffusion*grad2_omega
        flux_3 = -diffusion*grad3_omega


        call worker%integrate_volume_flux('Density * Omega','Diffusion',flux_1,flux_2,flux_3)

        !================================
        !       TURBULENCE FLUX - R_11
        !================================

        flux_1 = -diffusion_111
        flux_2 = -diffusion_112
        flux_3 = -diffusion_113


        call worker%integrate_volume_flux('Density * Reynolds-11','Diffusion',flux_1,flux_2,flux_3)

        !================================
        !       TURBULENCE FLUX - R_22
        !================================

        flux_1 = -diffusion_221
        flux_2 = -diffusion_222
        flux_3 = -diffusion_223


        call worker%integrate_volume_flux('Density * Reynolds-22','Diffusion',flux_1,flux_2,flux_3)

        !================================
        !       TURBULENCE FLUX - R_33
        !================================

        flux_1 = -diffusion_331
        flux_2 = -diffusion_332
        flux_3 = -diffusion_333


        call worker%integrate_volume_flux('Density * Reynolds-33','Diffusion',flux_1,flux_2,flux_3)


        !================================
        !       TURBULENCE FLUX - R_12
        !================================

        flux_1 = -diffusion_121
        flux_2 = -diffusion_122
        flux_3 = -diffusion_123


        call worker%integrate_volume_flux('Density * Reynolds-12','Diffusion',flux_1,flux_2,flux_3)
        !================================
        !       TURBULENCE FLUX - R_13
        !================================

        flux_1 = -diffusion_131
        flux_2 = -diffusion_132
        flux_3 = -diffusion_133


        call worker%integrate_volume_flux('Density * Reynolds-13','Diffusion',flux_1,flux_2,flux_3)

        !================================
        !       TURBULENCE FLUX - R_23
        !================================

        flux_1 = -diffusion_231
        flux_2 = -diffusion_232
        flux_3 = -diffusion_233


        call worker%integrate_volume_flux('Density * Reynolds-23','Diffusion',flux_1,flux_2,flux_3)


        !================================
        !       TURBULENCE FLUX - Energy 
        !================================

        flux_1 = -energy_flux_1
        flux_2 = -energy_flux_2
        flux_3 = -energy_flux_3


        call worker%integrate_volume_flux('Energy','Diffusion',flux_1,flux_2,flux_3)



    end subroutine compute
    !************************************************************************************************












end module rstm_ssglrrw_volume_diffusion
