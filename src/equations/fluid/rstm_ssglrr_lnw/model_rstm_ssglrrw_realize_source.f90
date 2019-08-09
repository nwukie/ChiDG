module model_rstm_ssglrrw_realize_source
#include <messenger.h>
    use mod_kinds,              only: rk
    use mod_constants,          only: ZERO, HALF, ONE, TWO, THREE, THIRD
    use type_model,             only: model_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    use mod_rstm_ssglrrw

    implicit none


    

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(model_t)  :: rstm_ssglrrw_realize_source_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_realize_source_t
    !***************************************************************************************





contains




    !>  Model to compute the realize_source term in the Reynolds Stress transport equation.
    !!  Note that this term is computed exactly from mean flow quantities with no modelling.
    !!
    !!  @author Eric M Wolf
    !!  @date   1/26/2018
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(rstm_ssglrrw_realize_source_t), intent(inout)   :: self

        call self%set_name('RSTMSSGLRRW Realizability Source')
        call self%set_dependency('f(Q-)')

        call self%add_model_field('Realizability Source-11')
        call self%add_model_field('Realizability Source-22')
        call self%add_model_field('Realizability Source-33')
        call self%add_model_field('Realizability Source-12')
        call self%add_model_field('Realizability Source-13')
        call self%add_model_field('Realizability Source-23')

        call self%add_model_field('Sustaining Source-R')
        call self%add_model_field('Sustaining Source-Omega')


    end subroutine init
    !***************************************************************************************






    !>  
    !!
    !!  @author Eric M Wolf
    !!  @date   1/26/2018
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(rstm_ssglrrw_realize_source_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable ::              &
                density, epsilon_t,                             &
                omega, R_sust, beta_w, omega_sust,              &
                r_11, r_22, r_33, r_12, r_13, r_23, det_r,      &
                src_11, src_22, src_33, src_12, src_13, src_23, &
                realizability_source_11, realizability_source_22, realizability_source_33, &
                realizability_source_12, realizability_source_13, realizability_source_23, tmp, &
                tmp_11, tmp_22, tmp_33, tmp_12, tmp_13, tmp_23
                
        real(rk), allocatable, dimension(:) :: weights, jinv

        density     = worker%get_field('Density',    'value')
        r_11 = worker%get_field('Density * Reynolds-11',    'value')/density
        r_22 = worker%get_field('Density * Reynolds-22',    'value')/density
        r_33 = worker%get_field('Density * Reynolds-33',    'value')/density
        r_12 = worker%get_field('Density * Reynolds-12',    'value')/density
        r_13 = worker%get_field('Density * Reynolds-13',    'value')/density
        r_23 = worker%get_field('Density * Reynolds-23',    'value')/density

        !det_R = density

        src_11 = ZERO*density
        src_22 = ZERO*density
        src_33 = ZERO*density
        src_12 = ZERO*density
        src_13 = ZERO*density
        src_23 = ZERO*density

        !src_11 = -SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_11*sin_ramp(-r_11, 0.0_rk, 0.01_rk*rstm_ssglrrw_k_infty)
        !src_22 = -SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_22*sin_ramp(-r_22, 0.0_rk, 0.01_rk*rstm_ssglrrw_k_infty)
        !src_33 = -SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_33*sin_ramp(-r_33, 0.0_rk, 0.01_rk*rstm_ssglrrw_k_infty)

        src_11 = -THIRD*SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_11*sin_ramp(-r_11, 0.0_rk, rstm_ssglrrw_R_infty)
        src_22 = -THIRD*SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_22*sin_ramp(-r_22, 0.0_rk, rstm_ssglrrw_R_infty)
        src_33 = -THIRD*SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_33*sin_ramp(-r_33, 0.0_rk, rstm_ssglrrw_R_infty)

        r_11 = worker%get_field('Reynolds-11',    'value')
        r_22 = worker%get_field('Reynolds-22',    'value')
        r_33 = worker%get_field('Reynolds-33',    'value')

        r_12 = worker%get_field('Reynolds-12',    'value')
        r_13 = worker%get_field('Reynolds-13',    'value')
        r_23 = worker%get_field('Reynolds-23',    'value')

        det_R = r_11*(r_22*r_33-r_23*r_23) - r_12*(r_12*r_33-r_23*r_13) + r_13*(r_12*r_23-r_22*r_13)

        src_12 = -THIRD*SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_12*(&
                    sin_ramp(r_12**TWO-(r_11*r_22), 0.0_rk, (rstm_ssglrrw_R_infty)**TWO) + &
                    sin_ramp(-det_R,0.0_rk, (rstm_ssglrrw_R_infty))**THREE)

        src_13 = -THIRD*SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_13*(&
                    sin_ramp(r_13**TWO-(r_11*r_33), 0.0_rk, (rstm_ssglrrw_R_infty)**TWO) + &
                    sin_ramp(-det_R,0.0_rk, (rstm_ssglrrw_R_infty)**THREE))

        src_23 = -THIRD*SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_23*(&
                    sin_ramp(r_23**TWO-(r_22*r_33), 0.0_rk, (rstm_ssglrrw_R_infty)**TWO) + &
                    sin_ramp(-det_R,0.0_rk, (rstm_ssglrrw_R_infty)**THREE))


        ! Try to averge the realizability source 
        !if (worker%interpolation_source == 'element') then
        !    weights = worker%quadrature_weights('element')
        !    jinv    = worker%inverse_jacobian('element')

        !    tmp_11 = src_11
        !    tmp_22 = src_22
        !    tmp_33 = src_33

        !    tmp_12 = src_12
        !    tmp_13 = src_13
        !    tmp_23 = src_23


        !    tmp_11 = sum(weights*jinv*src_11)/sum(weights*jinv)
        !    tmp_22 = sum(weights*jinv*src_22)/sum(weights*jinv)
        !    tmp_33 = sum(weights*jinv*src_33)/sum(weights*jinv)

        !    tmp_12 = sum(weights*jinv*src_12)/sum(weights*jinv)
        !    tmp_13 = sum(weights*jinv*src_13)/sum(weights*jinv)
        !    tmp_23 = sum(weights*jinv*src_23)/sum(weights*jinv)

        !    src_11 = tmp_11
        !    src_22 = tmp_22
        !    src_33 = tmp_33

        !    src_12 = tmp_12
        !    src_13 = tmp_13
        !    src_23 = tmp_23
        !end if


        !tmp = r_11
        !tmp = abs(r_11*r_22)
        !src_12 = -SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_12*&
        !            sin_ramp(r_12**TWO-tmp, 0.0_rk, (0.01_rk*rstm_ssglrrw_k_infty)**TWO) 

        !tmp = abs(r_11*r_33)
        !src_13 = -SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_13*&
        !            sin_ramp(r_13**TWO-tmp, 0.0_rk, (0.01_rk*rstm_ssglrrw_k_infty)**TWO) 

        !tmp = abs(r_22*r_33)
        !src_23 = -SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_23*&
        !            sin_ramp(r_23**TWO-tmp, 0.0_rk, (0.01_rk*rstm_ssglrrw_k_infty)**TWO) 

        !src_11 = -HALF*SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_11*sin_ramp(-r_11, 0.0_rk, 0.1_rk*rstm_ssglrrw_k_infty)
        !src_22 = -HALF*SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_22*sin_ramp(-r_22, 0.0_rk, 0.1_rk*rstm_ssglrrw_k_infty)
        !src_33 = -HALF*SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_33*sin_ramp(-r_33, 0.0_rk, 0.1_rk*rstm_ssglrrw_k_infty)


        !src_12 = -HALF*SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_12*(&
        !            sin_ramp(r_12**TWO-abs(r_11*r_22), 0.0_rk, (0.1_rk*rstm_ssglrrw_k_infty)**TWO) + &
        !            sin_ramp(-det_R,0.0_rk, (0.1_rk*rstm_ssglrrw_k_infty)**THREE))

        !src_13 = -HALF*SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_13*(&
        !            sin_ramp(r_13**TWO-abs(r_11*r_33), 0.0_rk, (0.1_rk*rstm_ssglrrw_k_infty)**TWO) + &
        !            sin_ramp(-det_R,0.0_rk, (0.1_rk*rstm_ssglrrw_k_infty)**THREE))

        !src_23 = -HALF*SSG_LRRW_cmu*rstm_ssglrrw_omega_infty*density*r_23*(&
        !            sin_ramp(r_23**TWO-abs(r_22*r_33), 0.0_rk, (0.1_rk*rstm_ssglrrw_k_infty)**TWO) + &
        !            sin_ramp(-det_R,0.0_rk, (0.1_rk*rstm_ssglrrw_k_infty)**THREE))



        call worker%store_model_field('Realizability Source-11', 'value', src_11)
        call worker%store_model_field('Realizability Source-22', 'value', src_22)
        call worker%store_model_field('Realizability Source-33', 'value', src_33)
        call worker%store_model_field('Realizability Source-12', 'value', src_12)
        call worker%store_model_field('Realizability Source-13', 'value', src_13)
        call worker%store_model_field('Realizability Source-23', 'value', src_23)


        R_sust = (TWO/THREE)*SSG_LRRW_cmu*density*rstm_ssglrrw_k_infty*rstm_ssglrrw_omega_infty
 

        omega     = worker%get_field('Density * Omega',    'value')/density
        beta_w    = worker%get_field('RSTMSSGLRRW Beta-w',    'value')
        omega_sust = beta_w*density*rstm_ssglrrw_omega_infty**TWO*exp(-omega)

        !if (worker%interpolation_source == 'element') then
        !    weights = worker%quadrature_weights('element')
        !    jinv    = worker%inverse_jacobian('element')

        !    tmp_11 = src_11
        !    tmp_22 = src_22

        !    tmp_11 = sum(weights*jinv*R_sust)/sum(weights*jinv)
        !    tmp_22 = sum(weights*jinv*omega_sust)/sum(weights*jinv)

        !    R_sust = tmp_11
        !    omega_sust  = tmp_22

        !end if


        call worker%store_model_field('Sustaining Source-R', 'value', R_sust)
        call worker%store_model_field('Sustaining Source-Omega', 'value', omega_sust)

    end subroutine compute
    !***************************************************************************************




end module model_rstm_ssglrrw_realize_source
