!>
!! Description: This model computes various quantities used in the SST TM.
!!
!! @author  Eric M. Wolf
!! @date    01/26/2018 
!!
!--------------------------------------------------------------------------------
module model_sst_turbulence_quantities
#include <messenger.h>
    use mod_kinds,              only: rk
    use mod_constants,          only: ZERO, HALF, ONE, TWO, THREE, FOUR
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
    type, extends(model_t)  :: sst_turbulence_quantities_t

    contains

        procedure   :: init
        procedure   :: compute

    end type sst_turbulence_quantities_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Eric M. Wolf
    !!  @date   1/26/2018
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(sst_turbulence_quantities_t), intent(inout)   :: self

        call self%set_name('SST Turbulence Quantities')
        call self%set_dependency('f(Grad(Q))')
        ! k
        ! gradient k
        call self%add_model_field('k - Gradient 1')
        call self%add_model_field('k - Gradient 2')
        call self%add_model_field('k - Gradient 3')



                ! omega
        call self%add_model_field('Omega')


        ! gradient omega
        call self%add_model_field('Omega - Gradient 1')
        call self%add_model_field('Omega - Gradient 2')
        call self%add_model_field('Omega - Gradient 3')
        call self%add_model_field('Omega Gradient Squared')



                !
        ! Turbulence Model Fields
        !
        call self%add_model_field('Turbulent Viscosity')
        call self%add_model_field('Second Coefficient of Turbulent Viscosity')
        call self%add_model_field('Turbulent Thermal Conductivity')

        ! SST source terms
        call self%add_model_field('SST Omega Source Term')
        call self%add_model_field('SST CD')
        call self%add_model_field('SST k Production Term')
        call self%add_model_field('SST k Destruction Term')
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
        class(sst_turbulence_quantities_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable ::              &
            density, mu,    &
            grad1_density,  grad2_density,  grad3_density,              &
            density_omega, omega,                                       &
            grad1_density_omega, grad2_density_omega, grad3_density_omega, &
            grad1_omega,    grad2_omega,    grad3_omega,                    &
            grad1_k, grad2_k, grad3_k,                            &
            mu_t, k, epsilon_t, temp1, temp2, omega_source_term, invdensity, grad_omega_sq, &
            alpha_star, beta_star, ReT, gmma, mu_l, &
            w_12, w_13, w_23, b_12, b_13, b_23, a_12, a_13, a_23, &
            tau_11, tau_22, tau_33, tau_12, tau_13, tau_23, &
            temp, omega_mod, src_term, cap_omega, &
            omega_lb, lamda_t, kap_t, &
            grad1_u, grad2_u, grad3_u, grad1_v, grad2_v, grad3_v, grad1_w, grad2_w, grad3_w, &
            grad1_density_k, grad2_density_k, grad3_density_k, &
            arg2, F2, distance, des_k, div_vel, &
            rot_12, rot_13, rot_23, str_11, str_22, str_33, str_12, str_13, str_23, chi_omega, f_beta, beta, divu, otemp1, otemp2

        real(rk), allocatable :: det_r(:)
        integer(ik) :: ii 


        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_field('Density',    'value')
        invdensity = ONE/density

        grad1_density = worker%get_field('Density', 'grad1')
        grad2_density = worker%get_field('Density', 'grad2')
        grad3_density = worker%get_field('Density', 'grad3')

        mu_l = worker%get_field('Laminar Viscosity',    'value')

    
        k = worker%get_field('Density * k',    'value')
        k = invdensity*k
        grad1_density_k = worker%get_field('Density * k', 'grad1')
        grad2_density_k = worker%get_field('Density * k', 'grad2')
        grad3_density_k = worker%get_field('Density * k', 'grad3')

        grad1_k = invdensity*(grad1_density_k-k*grad1_density)
        grad2_k = invdensity*(grad2_density_k-k*grad2_density)
        grad3_k = invdensity*(grad3_density_k-k*grad3_density)
        
        call worker%store_model_field('k - Gradient 1', 'value', grad1_k)
        call worker%store_model_field('k - Gradient 2', 'value', grad2_k)
        call worker%store_model_field('k - Gradient 3', 'value', grad3_k)

        omega = worker%get_field('Density * Omega',    'value')
        omega = omega*invdensity
        grad1_density_omega = worker%get_field('Density * Omega', 'grad1')
        grad2_density_omega = worker%get_field('Density * Omega', 'grad2')
        grad3_density_omega = worker%get_field('Density * Omega', 'grad3')


        grad1_omega = invdensity*(grad1_density_omega-omega*grad1_density)
        grad2_omega = invdensity*(grad2_density_omega-omega*grad2_density)
        grad3_omega = invdensity*(grad3_density_omega-omega*grad3_density)

        grad_omega_sq = grad1_omega*grad1_omega + grad2_omega*grad2_omega + grad3_omega*grad3_omega


        call worker%store_model_field('Omega - Gradient 1', 'value', grad1_omega)
        call worker%store_model_field('Omega - Gradient 2', 'value', grad2_omega)
        call worker%store_model_field('Omega - Gradient 3', 'value', grad3_omega)
        call worker%store_model_field('Omega Gradient Squared', 'value', grad_omega_sq)
 



        !
        ! Compute lower bound for Omega to enforce realizability
        !

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


        a_12 = str_11+str_22
        a_13 = str_11+str_33
        a_23 = str_22+str_33
        !print *, 'a_12'
        !print *, a_12(:)%x_ad_
        !print *, 'a_13'
        !print *, a_13(:)%x_ad_
        !print *, 'a_23'
        !print *, a_23(:)%x_ad_


        b_12 = str_11*str_22-str_12**TWO
        b_13 = str_11*str_33-str_13**TWO
        b_23 = str_22*str_33-str_23**TWO
        !print *, 'b_12'
        !print *, b_12(:)%x_ad_
        !print *, 'b_13'
        !print *, b_13(:)%x_ad_
        !print *, 'b_23'
        !print *, b_23(:)%x_ad_

        w_12 = a_12
        w_13 = a_13
        w_23 = a_23

        w_12 = 1.5_rk*(0.5_rk*a_12+sqrt(abs((0.5_rk*a_12)**TWO-b_12) + 1.0e-12_rk))
        w_13 = 1.5_rk*(0.5_rk*a_13+sqrt(abs((0.5_rk*a_13)**TWO-b_13) + 1.0e-12_rk))
        w_23 = 1.5_rk*(0.5_rk*a_23+sqrt(abs((0.5_rk*a_23)**TWO-b_23) + 1.0e-12_rk))

        !print *, 'w_12'
        !print *, w_12(:)%x_ad_
        !print *, 'w_13'
        !print *, w_13(:)%x_ad_
        !print *, 'w_23'
        !print *, w_23(:)%x_ad_

        temp = ZERO*density
        omega_lb = ZERO*density
        temp = max(w_12, w_13)
        
        omega_lb = max(temp, w_23)
        !print *, 'omega_lb'
        !print *, omega_lb(:)%x_ad_


        ! NOTE: imposing the lower bound on omega is resulting in 
        ! NaN residual in first Quasi-Newton step???
        omega_mod = omega
        !omega_mod = log(max(exp(omega), omega_lb))
        otemp1 = omega
        otemp2 = omega
        otemp1 = exp(omega)
        otemp2 = max(otemp1, omega_lb)
        !print *, 'otemp2'
        !print *, otemp2(:)%x_ad_
        omega_mod = log(otemp2)
        !omega_mod(:)%xp_ad_(:) = omega(:)%xp_ad_(:)
        !print *, 'omega mod'
        !print *, omega_mod(:)%x_ad_
        !omega_mod = omega

        call worker%store_model_field('Omega', 'value', omega_mod)

                
        omega_source_term = density*exp(-omega_mod)*(grad1_k*grad1_omega + grad2_k*grad2_omega + grad3_k*grad3_omega)

        call worker%store_model_field('SST Omega Source Term', 'value', omega_source_term)

        omega_source_term = density*(grad1_k*grad1_omega + grad2_k*grad2_omega + grad3_k*grad3_omega)

        call worker%store_model_field('SST CD', 'value', omega_source_term)


        k = worker%get_field('k',    'value')
        if (worker%interpolation_source == 'boundary') then
            distance            = worker%get_field('Wall Distance',             'value', 'face interior') 
        else if (worker%interpolation_source == 'face exterior') then
            distance            = worker%get_field('Wall Distance',             'value', 'face interior') 
        else
            distance            = worker%get_field('Wall Distance',             'value') 
        end if


        

        ! Compute blending factor
        arg2 = max(&
                    TWO*sqrt(k + 1.0e-12_rk)/(0.09_rk*exp(omega_mod)*abs(distance)+1.0e-12_rk), 500.0_rk*mu_l/(density*exp(omega_mod)*abs(distance)**TWO+1.0e-12_rk)&
        )
        F2 = arg2
        arg2 = F2**TWO
        temp1 = F2
        temp1 = ONE-exp(-TWO*arg2)
        temp2 = temp1
        temp2 = exp(-TWO*arg2)+ONE
        F2 = temp1/temp2
        !F2 = tanh(arg2(:)%x_ad_**TWO)
        !print *, 'F2'
        !print *, F2(:)%x_ad_
       
        rot_12 = worker%get_field('Rotation Rate-12', 'value')
        rot_13 = worker%get_field('Rotation Rate-13', 'value')
        rot_23 = worker%get_field('Rotation Rate-23', 'value')

        cap_omega = omega
        cap_omega = TWO*sqrt(rot_12**TWO+rot_13**TWO+rot_23**TWO + 1.0e-12_rk)


        temp = max(0.31_rk*exp(omega_mod), F2*cap_omega)

        mu_t = 0.31_rk*density*k/(temp + 1.0e-12_rk)
        !mu_t = density*k*exp(-omega_mod)
        lamda_t = (-TWO/THREE)*mu_t
        kap_t = cp*mu_t/0.9_rk


        call worker%store_model_field('Turbulent Viscosity',                       'value', mu_t   )
        call worker%store_model_field('Second Coefficient of Turbulent Viscosity', 'value', lamda_t)
        call worker%store_model_field('Turbulent Thermal Conductivity',            'value', kap_t    )


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

        src_term =  tau_11*grad1_u + tau_12*grad2_u + tau_13*grad3_u + &
                    tau_12*grad1_v + tau_22*grad2_v + tau_23*grad3_v + &
                    tau_13*grad1_w + tau_23*grad2_w + tau_33*grad3_w

        call worker%store_model_field('SST k Production Term', 'value', src_term)


        des_k = density*k*exp(omega_mod)
    
        call worker%store_model_field('SST k Destruction Term', 'value', des_k)

    end subroutine compute
    !***************************************************************************************




end module model_sst_turbulence_quantities
