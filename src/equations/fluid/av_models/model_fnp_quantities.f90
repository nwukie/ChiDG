module model_fnp_quantities
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: HALF, ONE, TWO, PI, RKTOL
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    use ieee_arithmetic 
    implicit none


    


    !> Based on the reference:
    !! Fernandez, Pablo, Cuong Nguyen, and Jaime Peraire. 
    !! "A physics-based shock capturing method for unsteady laminar and turbulent flows." 
    !! 2018 AIAA Aerospace Sciences Meeting. 2018.
    !!
    !! @author  Eric M. Wolf
    !! @date    09/06/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(model_t)  :: fnp_quantities_t


    contains

        procedure   :: init
        procedure   :: compute

    end type fnp_quantities_t
    !***************************************************************************************





contains




    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    09/06/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)   
        class(fnp_quantities_t), intent(inout)   :: self

        call self%set_name('FNP Quantities')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Sound Speed at Critical Temperature')
        call self%add_model_field('Temperature Gradient Under Reference Metric')
        call self%add_model_field('Stagnation Temperature')
        call self%add_model_field('Spectral Norm Deleted Diagonal Velocity Gradient')
        call self%add_model_field('Maximum Isentropic Velocity')
    end subroutine init
    !***************************************************************************************






    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    07/11/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker)
        class(fnp_quantities_t),  intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        integer(ik) :: ii
        type(AD_D), allocatable,    dimension(:) ::         &
            grad1_u,       grad2_u,       grad3_u,          &
            grad1_v,       grad2_v,       grad3_v,          &
            grad1_w,       grad2_w,       grad3_w,          &
            curl_vel_1, curl_vel_2, curl_vel_3,             &
            s_omega, s_d,                                   &
            s_beta, s_beta_smoothed,                        &
            s_kappa, s_kappa_smoothed,                        &
            s_mu, s_mu_smoothed,                        &
            s_thr, s_min, s_max                        &
            t_xi, t_stag, c_star, Lv_norm, v_max,   &
            div_vel, density, mom1, mom2, mom3, vel_mag, cutoff, D, Dm

        real(rk),   allocatable,    dimension(:) :: h
        real(rk) :: Dmax, d0, delta_d

        density = worker%get_field('Density',    'value')
        mom1    = worker%get_field('Momentum-1', 'value')
        mom2    = worker%get_field('Momentum-2', 'value')
        mom3    = worker%get_field('Momentum-3', 'value')

        ! The reference describes "sound speed at critical temperature"
        ! Critical temperature of air is given as a constant T_star = 132.63 K.
        ! But in the reference T_star is said to have spatial variation?
        c_star = sqrt(gam*R*132.63_rk)

        call worker%store_model_field('Sound Speed at Critical Temperature', 'value', c_star)

        temp_grad1 = worker%get_field('Temperature Gradient - 1', 'value')
        temp_grad2 = worker%get_field('Temperature Gradient - 2', 'value')
        temp_grad3 = worker%get_field('Temperature Gradient - 3', 'value')

        ! Note:
        dxidx = worker%get_metric('deformed')
        do inode = 1, nnodes
            metric(:,:,inode) = inv_3x3(dxidx(:,:,inode))

        end do

        temp_grad_mod1 = temp_grad1*metric(1,1,:) + temp_grad2*metric(1,2,:) +  temp_grad3*metric(1,3,:)
        temp_grad_mod2 = temp_grad1*metric(2,1,:) + temp_grad2*metric(2,2,:) +  temp_grad3*metric(2,3,:)
        temp_grad_mod3 = temp_grad1*metric(3,1,:) + temp_grad2*metric(3,2,:) +  temp_grad3*metric(3,3,:)

        temp_grad_mod_norm = temp_grad_mod_1
        temp_grad_mod_norm = sqrt(temp_grad_mod1**TWO+ temp_grad_mod2**TWO + temp_grad_mod3**TWO)
        call worker%store_model_field('Temperature Gradient Under Reference Metric - Magnitude', 'value', temp_grad_mod_norm) 

        t = worker%get_field('Temperature', 'value') 

        vel_mag = (mom1**TWO+mom2**TWO+mom3**TWO)/(density**TWO )

        t_stag = t + vel_mag/(TWO*cp)

        call worker%store_model_field('Stagnation Temperature', 'value', t_stag)


        vel1_grad1 = worker%get_field('Velocity Gradient 1 - 1', 'value')
        vel1_grad2 = worker%get_field('Velocity Gradient 1 - 2', 'value')
        vel1_grad3 = worker%get_field('Velocity Gradient 1 - 3', 'value')

        vel2_grad1 = worker%get_field('Velocity Gradient 2 - 1', 'value')
        vel2_grad2 = worker%get_field('Velocity Gradient 2 - 2', 'value')
        vel2_grad3 = worker%get_field('Velocity Gradient 2 - 3', 'value')

        vel3_grad1 = worker%get_field('Velocity Gradient 3 - 1', 'value')
        vel3_grad2 = worker%get_field('Velocity Gradient 3 - 2', 'value')
        vel3_grad3 = worker%get_field('Velocity Gradient 3 - 3', 'value')


        !Instead of the spectral norm (2-norm), use the Frobenius norm for ease of computation
        do inode = 1, nnodes
            frob_norm(inode) = &
            vel1_grad2(inode)**TWO + vel1_grad3(inode)**TWO + &
            vel2_grad1(inode)**TWO + vel2_grad3(inode)**TWO + &
            vel3_grad2(inode)**TWO + vel3_grad1(inode)**TWO 

        end do

        frob_norm = sqrt(frob_norm)

        call worker%store_model_field('Spectral Norm Deleted Diagonal Velocity Gradient', 'value', frob_norm)

        c = worker%get_field('Sound Speed', 'value')

        vel_max = c

        vel_max = sqrt(vel_mag + (TWO/(gam-ONE))*c**TWO)

        call worker%store_model_field('Maximum Isentropic Velocity', 'value', vel_max)



        density_grad1 = worker%get_field('Density', 'grad1')
        density_grad2 = worker%get_field('Density', 'grad2')
        density_grad3 = worker%get_field('Density', 'grad3')

        do inode = 1, nnodes

            dmd(inode) = &
            density_grad1(inode)*(&
            metric(1,1)*density_grad1(inode) + metric(1,2)*density_grad2(inode) + metric(1,3)*density_grad3(inode)
            )+&
            density_grad2(inode)*(&
            metric(2,1)*density_grad1(inode) + metric(2,2)*density_grad2(inode) + metric(2,3)*density_grad3(inode)
            )+&
            density_grad3(inode)*(&
            metric(3,1)*density_grad1(inode) + metric(3,2)*density_grad2(inode) + metric(3,3)*density_grad3(inode)
            )

            tmt(inode) = &
            temp_grad1(inode)*(&
            metric(1,1)*temp_grad1(inode) + metric(1,2)*temp_grad2(inode) + metric(1,3)*temp_grad3(inode)
            )+&
            temp_grad2(inode)*(&
            metric(2,1)*temp_grad1(inode) + metric(2,2)*temp_grad2(inode) + metric(2,3)*temp_grad3(inode)
            )+&
            temp_grad3(inode)*(&
            metric(3,1)*temp_grad1(inode) + metric(3,2)*temp_grad2(inode) + metric(3,3)*temp_grad3(inode)
            )

        end do

        h_beta = h_ref*(sqrt(density_grad1**TWO + density_grad2**TWO + density_grad3**TWO )/&
        sqrt(dmd + RKTOL**TWO ))


        h_kappa = h_ref*(sqrt(temp_grad1**TWO + temp_grad2**TWO + temp_grad3**TWO )/&
        sqrt(tmt + RKTOL**TWO ))
    end subroutine compute
    !***************************************************************************************




end module model_fnp_quantities
