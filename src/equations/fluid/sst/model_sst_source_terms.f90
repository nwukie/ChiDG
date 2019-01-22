!>
!! Description: This model computes source terms for the SST TM.
!!
!! @author  Eric M. Wolf
!! @date    01/26/2018 
!!
!--------------------------------------------------------------------------------
module model_sst_source_terms
#include <messenger.h>
    use mod_kinds,              only: rk
    use mod_constants,          only: ZERO, HALF, ONE, TWO, THREE, FOUR, PI
    use type_model,             only: model_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    use ieee_arithmetic,        only: ieee_is_nan
    use mod_sst
    use mod_fluid

    implicit none


    

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(model_t)  :: sst_source_terms_t

        real(rk)    :: k_infty      = 1.0e-3_rk
        real(rk)    :: omega_infty  = 345.0_rk

    contains

        procedure   :: init
        procedure   :: compute

    end type sst_source_terms_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Eric M. Wolf
    !!  @date   1/26/2018
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(sst_source_terms_t), intent(inout)   :: self
        real(rk)            :: k_infty, omega_infty
        integer             :: unit, msg
        logical             :: file_exists

        namelist /k_omega/   k_infty, omega_infty 




        call self%set_name('SST Source Terms')
        call self%set_dependency('f(Grad(Q))')

        ! SST source terms
        call self%add_model_field('SST Omega Source Term')
        call self%add_model_field('SST k Source Term')
        call self%add_model_field('SST Energy Source Term')
!
        ! Check if input from 'models.nml' is available.
        !   1: if available, read and set self%k_infty/omega_infty
        !   2: if not available, do nothing and retain default values
        !
        inquire(file='models.nml', exist=file_exists)
        if (file_exists) then
            open(newunit=unit,form='formatted',file='models.nml')
            read(unit,nml=k_omega,iostat=msg)
            if (msg == 0) self%k_infty = k_infty 
            if (msg == 0) self%omega_infty = omega_infty 
            close(unit)
        end if


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
        class(sst_source_terms_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable ::              &
            density, mu,    &
            alpha_w, beta_w, beta_k, sigma_w, sigma_d, &
            p_ub, production, production_mod, destruction, destruction_neg, &
            CD, k_src, omega_src, omega_prod_fac,&
            k_sust, omega_sust, &
            grad1_density,  grad2_density,  grad3_density,              &
            density_omega, omega,                                       &
            grad1_density_omega, grad2_density_omega, grad3_density_omega, &
            grad1_omega,    grad2_omega,    grad3_omega,                    &
            grad1_k, grad2_k, grad3_k,                            &
            mu_t, density_k, k, epsilon_t, temp1, temp2, omega_source_term, invdensity, grad_omega_sq, &
            alpha_star, beta_star, ReT, gmma, mu_l, &
            w_12, w_13, w_23, b_12, b_13, b_23, a_12, a_13, a_23, &
            tau_11, tau_22, tau_33, tau_12, tau_13, tau_23, &
            temp, omega_mod, src_term, cap_omega, &
            omega_lb, lamda_t, kap_t, &
            grad1_u, grad2_u, grad3_u, grad1_v, grad2_v, grad3_v, grad1_w, grad2_w, grad3_w, &
            grad1_density_k, grad2_density_k, grad3_density_k, &
            arg2, F2, distance, des_k, div_vel, e_src,&
            rot_12, rot_13, rot_23, str_11, str_22, str_33, str_12, str_13, str_23, chi_omega, f_beta, beta, divu, otemp1, otemp2

        real(rk), allocatable :: det_r(:)
        real(rk), allocatable :: dummy(:)
        real(rk) :: b, c
        integer(ik) :: ii 


        !
        ! Interpolate solution to quadrature nodes
        !
        density     = worker%get_field('Density',    'value')
        density_k     = worker%get_field('Density * k',    'value')
        invdensity  = ONE/density

        mu_l = worker%get_field('Laminar Viscosity',    'value')
        mu_t = worker%get_field('Turbulent Viscosity',    'value')

    
        ! SST coefficients
        alpha_w = worker%get_field('SST alpha_w', 'value')
        beta_w  = worker%get_field('SST beta_w', 'value')
        beta_k  = worker%get_field('SST beta_k', 'value')
        sigma_w = worker%get_field('SST sigma_w', 'value')
        sigma_d = worker%get_field('SST sigma_d', 'value')

        k       = worker%get_field('k',    'value')
        omega   = worker%get_field('Omega',    'value')

        str_11 = TWO*worker%get_field('Strain Rate-11', 'value')
        str_22 = TWO*worker%get_field('Strain Rate-22', 'value')
        str_33 = TWO*worker%get_field('Strain Rate-33', 'value')

        str_12 = TWO*worker%get_field('Strain Rate-12', 'value')
        str_13 = TWO*worker%get_field('Strain Rate-13', 'value')
        str_23 = TWO*worker%get_field('Strain Rate-23', 'value')
        
        div_vel = str_11
        div_vel = HALF*(str_11 + str_22 + str_33)
        str_11 = str_11 - (TWO/THREE)*div_vel
        str_22 = str_22 - (TWO/THREE)*div_vel
        str_33 = str_33 - (TWO/THREE)*div_vel


        grad1_u = worker%get_field('Velocity 1 - Gradient 1', 'value')
        grad2_u = worker%get_field('Velocity 1 - Gradient 2', 'value')
        grad3_u = worker%get_field('Velocity 1 - Gradient 3', 'value')

        grad1_v = worker%get_field('Velocity 2 - Gradient 1', 'value')
        grad2_v = worker%get_field('Velocity 2 - Gradient 2', 'value')
        grad3_v = worker%get_field('Velocity 2 - Gradient 3', 'value')
        
        grad1_w = worker%get_field('Velocity 3 - Gradient 1', 'value')
        grad2_w = worker%get_field('Velocity 3 - Gradient 2', 'value')
        grad3_w = worker%get_field('Velocity 3 - Gradient 3', 'value')

        tau_11 = mu_t*str_11 - (TWO/THREE)*density*k
        tau_22 = mu_t*str_22 - (TWO/THREE)*density*k
        tau_33 = mu_t*str_33 - (TWO/THREE)*density*k

        tau_12 = mu_t*str_12
        tau_13 = mu_t*str_13
        tau_23 = mu_t*str_23

        production =    tau_11*grad1_u + tau_12*grad2_u + tau_13*grad3_u + &
                        tau_12*grad1_v + tau_22*grad2_v + tau_23*grad3_v + &
                        tau_13*grad1_w + tau_23*grad2_w + tau_33*grad3_w

        destruction = beta_k*density*k*exp(omega)

        ! Extra destruction term in negative k regions
        allocate(dummy(size(destruction)))
        dummy = ZERO
        !destruction_neg = beta_k*sst_omega_infty*density_k*sin_ramp(-invdensity*density_k,0.0_rk*dummy,10.0_rk*sst_k_infty + 0.0_rk*dummy)/sst_omega_infty
        destruction_neg = beta_k*sst_omega_infty*density_k*sin_ramp(-invdensity*density_k,0.0_rk,10.0_rk*sst_k_infty)/sst_omega_infty
        !destruction_neg = ZERO

        ! Limit production term
        ! Use a soft min instead of hard min:
        !production_mod = min(production, 10.0_rk*destruction)
        
        production_mod = production
        !p_ub = 10.0_rk*destruction

        !b = 100.0_rk
        !c = HALF - atan(b)/PI
        !production_mod = p_ub - (p_ub - production)*(atan(b*(p_ub-production))/PI + HALF) - c
        !production_mod = smax(p_ub, production)
        
        k_sust = beta_k*density*self%k_infty*self%omega_infty
        k_src = production_mod - destruction - destruction_neg + k_sust

        call worker%store_model_field('SST k Source Term', 'value', k_src)

        e_src = -production_mod + destruction - k_sust
        call worker%store_model_field('SST Energy Source Term', 'value', e_src)

        ! Compute the omega production term in an algebraically simplified form
        omega_prod_fac =  worker%get_field('SST Omega Production Simplified Factor', 'value')
        tau_11 = str_11 - (TWO/THREE)*omega_prod_fac
        tau_22 = str_22 - (TWO/THREE)*omega_prod_fac
        tau_33 = str_33 - (TWO/THREE)*omega_prod_fac

        tau_12 = str_12
        tau_13 = str_13
        tau_23 = str_23

        production =    tau_11*grad1_u + tau_12*grad2_u + tau_13*grad3_u + &
                        tau_12*grad1_v + tau_22*grad2_v + tau_23*grad3_v + &
                        tau_13*grad1_w + tau_23*grad2_w + tau_33*grad3_w



        grad_omega_sq =  worker%get_field('Omega Gradient Squared', 'value')
        CD =  worker%get_field('SST CD', 'value')
        omega_sust = beta_w*density*self%omega_infty**TWO*exp(-omega)
        omega_src = (alpha_w*density*exp(-omega))*production &
                        - beta_w*density*exp(omega) &
                        + sigma_d*CD &
                        + (mu_l + sigma_w*mu_t)*grad_omega_sq  + omega_sust



        
    
        call worker%store_model_field('SST Omega Source Term', 'value', omega_src)


    end subroutine compute
    !***************************************************************************************




end module model_sst_source_terms
