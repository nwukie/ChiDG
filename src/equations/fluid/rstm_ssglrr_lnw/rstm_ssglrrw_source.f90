module rstm_ssglrrw_source
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,SIX,HALF,PI
    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    implicit none

    private

    
    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/31/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: rstm_ssglrrw_source_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_source_operator_t
    !******************************************************************************










contains

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/31/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rstm_ssglrrw_source_operator_t),   intent(inout)      :: self

        ! Set operator name.
        call self%set_name('RSTMSSGLRRW Source Operator')

        ! Set operator type.
        call self%set_operator_type('Volume Diffusive Operator')

        ! Set operator equations being integrated.
        call self%add_primary_field('Density * Omega')
        call self%add_primary_field('Density * Reynolds-11')
        call self%add_primary_field('Density * Reynolds-22')
        call self%add_primary_field('Density * Reynolds-33')
        call self%add_primary_field('Density * Reynolds-12')
        call self%add_primary_field('Density * Reynolds-13')
        call self%add_primary_field('Density * Reynolds-23')

        call self%add_auxiliary_field('Wall Distance : p-Poisson')
        call self%add_model('Wall Distance : p-Poisson Normalization')
    end subroutine init
    !********************************************************************************


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/31/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(rstm_ssglrrw_source_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:) ::                            &
            density, density_omega, production_trace, omega, k_t,           &
            alpha_w, beta_w, sigma_d_w, omega_source_term,                                     &
            production_11, production_22, production_33,                    &
            production_12, production_13, production_23,                    &
            pressure_strain_11, pressure_strain_22, pressure_strain_33,                    &
            pressure_strain_12, pressure_strain_13, pressure_strain_23,                    &
            dissipation_11, dissipation_22, dissipation_33,                    &
            dissipation_12, dissipation_13, dissipation_23,                    &
            rsrc_11, rsrc_22, rsrc_33, rsrc_12, rsrc_13, rsrc_23, &
            source, mu_l, mu_t

        real(rk)    :: const, epsilon_vorticity, eps

        ! Omega source = production_trace_term - omega_squared_term + omega_source_term

        density_omega       = worker%get_field('Density * Omega',           'value', 'element')
        production_trace    = worker%get_field('Production-Trace',          'value', 'element')
        omega               = worker%get_field('Omega',                     'value', 'element')
        omega_source_term   = worker%get_field('Omega Gradient Squared',                     'value', 'element')
        k_t                 = worker%get_field('Turbulence Kinetic Energy', 'value', 'element')
        alpha_w             = worker%get_field('RSTMSSGLRRW Alpha-w',       'value', 'element')
        beta_w              = worker%get_field('RSTMSSGLRRW Beta-w',        'value', 'element')
        sigma_d_w           = worker%get_field('RSTMSSGLRRW Sigma_d-w',     'value', 'element')

        mu_l       = worker%get_field('Laminar Viscosity',           'value', 'element')
        mu_t       = worker%get_field('Equivalent Eddy Viscosity',           'value', 'element')
        ! Reynold stress source = production+pressure_strain-dissipation

        production_11   = worker%get_field('Production-11', 'value', 'element')
        production_22   = worker%get_field('Production-22', 'value', 'element')
        production_33   = worker%get_field('Production-33', 'value', 'element')
        production_12   = worker%get_field('Production-12', 'value', 'element')
        production_13   = worker%get_field('Production-13', 'value', 'element')
        production_23   = worker%get_field('Production-23', 'value', 'element')

        pressure_strain_11  = worker%get_field('Pressure-Strain-11', 'value', 'element')
        pressure_strain_22  = worker%get_field('Pressure-Strain-22', 'value', 'element')
        pressure_strain_33  = worker%get_field('Pressure-Strain-33', 'value', 'element')
        pressure_strain_12  = worker%get_field('Pressure-Strain-12', 'value', 'element')
        pressure_strain_13  = worker%get_field('Pressure-Strain-13', 'value', 'element')
        pressure_strain_23  = worker%get_field('Pressure-Strain-23', 'value', 'element')

        dissipation_11      = worker%get_field('Turbulence Dissipation-11', 'value', 'element')
        dissipation_22      = worker%get_field('Turbulence Dissipation-22', 'value', 'element')
        dissipation_33      = worker%get_field('Turbulence Dissipation-33', 'value', 'element')
        dissipation_12      = worker%get_field('Turbulence Dissipation-12', 'value', 'element')
        dissipation_13      = worker%get_field('Turbulence Dissipation-13', 'value', 'element')
        dissipation_23      = worker%get_field('Turbulence Dissipation-23', 'value', 'element')

        rsrc_11   = worker%get_field('Realizability Source-11', 'value', 'element')
        rsrc_22   = worker%get_field('Realizability Source-22', 'value', 'element')
        rsrc_33   = worker%get_field('Realizability Source-33', 'value', 'element')
        rsrc_12   = worker%get_field('Realizability Source-12', 'value', 'element')
        rsrc_13   = worker%get_field('Realizability Source-13', 'value', 'element')
        rsrc_23   = worker%get_field('Realizability Source-23', 'value', 'element')
        !
        ! Interpolate solution to quadrature nodes
        !
        density     = worker%get_field('Density'          , 'value', 'element')



        !========================================================================
        !                       Omega Source Term
        !========================================================================
        !source = HALF*alpha_w*omega*production_trace/k_t - beta_w*density_omega*omega+sigma_d_w*omega_source_term
        source = HALF*(5.0_rk/9.0_rk)*production_trace/(k_t+1.0e-16) -(3.0_rk/40.0_rk)*density*exp(omega) + (mu_l+0.5_rk*mu_t)*omega_source_term
        !source = ZERO

        call worker%integrate_volume_source('Density * Omega',source)

        !========================================================================
        !                       R_11 Source Term
        !========================================================================
        source = production_11 + pressure_strain_11 - dissipation_11 + rsrc_11
        !source = ZERO

        call worker%integrate_volume_source('Density * Reynolds-11',source)

        !========================================================================
        !                       R_22 Source Term
        !========================================================================
        source = production_22 + pressure_strain_22 - dissipation_22 + rsrc_22
        !source = ZERO

        call worker%integrate_volume_source('Density * Reynolds-22',source)

        !========================================================================
        !                       R_33 Source Term
        !========================================================================
        source = production_33 + pressure_strain_33 - dissipation_33 + rsrc_33
        !source = ZERO

        call worker%integrate_volume_source('Density * Reynolds-33',source)

        !========================================================================
        !                       R_12 Source Term
        !========================================================================
        source = production_12 + pressure_strain_12 - dissipation_12 + rsrc_12
        !source = ZERO

        call worker%integrate_volume_source('Density * Reynolds-12',source)

        !========================================================================
        !                       R_13 Source Term
        !========================================================================
        source = production_13 + pressure_strain_13 - dissipation_13 + rsrc_13
        !source = ZERO

        call worker%integrate_volume_source('Density * Reynolds-13',source)

        !========================================================================
        !                       R_23 Source Term
        !========================================================================
        source = production_23 + pressure_strain_23 - dissipation_23 + rsrc_23
        !source = ZERO

        call worker%integrate_volume_source('Density * Reynolds-23',source)


    end subroutine compute
    !*********************************************************************************************************






end module rstm_ssglrrw_source
