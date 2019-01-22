module model_rstm_ssglrrw_artificial_viscosity
#include <messenger.h>
    use mod_kinds,              only: rk
    use mod_constants,          only: ZERO, HALF, ONE, TWO, THREE
    use mod_rstm_ssglrrw
    use type_model,             only: model_t
    use type_chidg_worker,      only: chidg_worker_t
    use ieee_arithmetic,        only: ieee_is_nan
    use DNAD_D

    implicit none


    

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(model_t)  :: rstm_ssglrrw_artificial_viscosity_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_artificial_viscosity_t
    !***************************************************************************************





contains



    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)   
        class(rstm_ssglrrw_artificial_viscosity_t), intent(inout)   :: self

        call self%set_name('RSTMSSGLRRW Artificial Viscosity')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('RSTM AV-11')
        call self%add_model_field('RSTM AV-22')
        call self%add_model_field('RSTM AV-33')


    end subroutine init
    !***************************************************************************************





    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(rstm_ssglrrw_artificial_viscosity_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable ::              &
            density_reynolds_11, density_reynolds_22, density_reynolds_33, &
            density_reynolds_12, density_reynolds_13, density_reynolds_23, &
            mu_neg_11, mu_neg_22, mu_neg_33, &
            invdensity, density               

        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_field('Density',    'value')
        invdensity = ONE/density

        density_reynolds_11 = worker%get_field('Density * Reynolds-11',    'value')
        density_reynolds_22 = worker%get_field('Density * Reynolds-22',    'value')
        density_reynolds_33 = worker%get_field('Density * Reynolds-33',    'value')
        density_reynolds_12 = worker%get_field('Density * Reynolds-12',    'value')
        density_reynolds_13 = worker%get_field('Density * Reynolds-13',    'value')
        density_reynolds_23 = worker%get_field('Density * Reynolds-23',    'value')


        mu_neg_11 = -density_reynolds_11*sin_ramp(-invdensity*density_reynolds_11, 0.0_rk, 10.0_rk*rstm_ssglrrw_k_infty)/rstm_ssglrrw_omega_infty
        mu_neg_22 = -density_reynolds_22*sin_ramp(-invdensity*density_reynolds_22, 0.0_rk, 10.0_rk*rstm_ssglrrw_k_infty)/rstm_ssglrrw_omega_infty
        mu_neg_33 = -density_reynolds_33*sin_ramp(-invdensity*density_reynolds_33, 0.0_rk, 10.0_rk*rstm_ssglrrw_k_infty)/rstm_ssglrrw_omega_infty

        call worker%store_model_field('RSTM AV-11','value',mu_neg_11) 
        call worker%store_model_field('RSTM AV-22','value',mu_neg_22) 
        call worker%store_model_field('RSTM AV-33','value',mu_neg_33) 
    end subroutine compute
    !***************************************************************************************




end module model_rstm_ssglrrw_artificial_viscosity
