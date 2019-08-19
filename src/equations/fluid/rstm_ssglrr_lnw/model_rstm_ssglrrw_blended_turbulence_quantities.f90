!>
!! Description: This model computes various quantities used in the SSG-LRR-w RSM.
!!
!! @author  Eric M. Wolf
!! @date    01/26/2018 
!!
!--------------------------------------------------------------------------------
module model_rstm_ssglrrw_blended_turbulence_quantities
#include <messenger.h>
    use mod_kinds,              only: rk
    use mod_constants,          only: ZERO, HALF, ONE, TWO, THREE
    use type_model,             only: model_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    use ieee_arithmetic,        only: ieee_is_nan
    use mod_fluid,              only: cp
    use mod_rstm_ssglrrw

    implicit none


    

    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    01/26/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(model_t)  :: rstm_ssglrrw_blended_turbulence_quantities_t

        real(rk)    :: blend = 0.0_rk
    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_blended_turbulence_quantities_t
    !***************************************************************************************





contains




    !>  Initialize the model with a name and the model fields it is contributing to.
    !!
    !!  @author Eric M. Wolf
    !!  @date   1/26/2018
    !!
    !---------------------------------------------------------------------------------------
    subroutine init(self)   
        class(rstm_ssglrrw_blended_turbulence_quantities_t), intent(inout)   :: self
        real(rk)            :: blend 
        integer             :: unit, msg
        logical             :: file_exists

        namelist /blended_rstm/   blend 


        call self%set_name('RSTMSSGLRRW Blended Turbulence Quantities')
        call self%set_dependency('f(Grad(Q))')

                ! omega
        call self%add_model_field('Omega')

        ! gradient omega
        call self%add_model_field('Omega - Gradient 1')
        call self%add_model_field('Omega - Gradient 2')
        call self%add_model_field('Omega - Gradient 3')
        call self%add_model_field('Omega Gradient Squared')


        !call self%add_model_field('Reynolds-11')
        !call self%add_model_field('Reynolds-22')
        !call self%add_model_field('Reynolds-33')
        !call self%add_model_field('Reynolds-12')
        !call self%add_model_field('Reynolds-13')
        !call self%add_model_field('Reynolds-23')

        !call self%add_model_field('Reynolds-Stress-11')
        !call self%add_model_field('Reynolds-Stress-22')
        !call self%add_model_field('Reynolds-Stress-33')
        !call self%add_model_field('Reynolds-Stress-12')
        !call self%add_model_field('Reynolds-Stress-13')
        !call self%add_model_field('Reynolds-Stress-23')




        call self%add_model_field('Reynolds-11 - Gradient 1')
        call self%add_model_field('Reynolds-11 - Gradient 2')
        call self%add_model_field('Reynolds-11 - Gradient 3')

        call self%add_model_field('Reynolds-22 - Gradient 1')
        call self%add_model_field('Reynolds-22 - Gradient 2')
        call self%add_model_field('Reynolds-22 - Gradient 3')

        call self%add_model_field('Reynolds-33 - Gradient 1')
        call self%add_model_field('Reynolds-33 - Gradient 2')
        call self%add_model_field('Reynolds-33 - Gradient 3')

        call self%add_model_field('Reynolds-12 - Gradient 1')
        call self%add_model_field('Reynolds-12 - Gradient 2')
        call self%add_model_field('Reynolds-12 - Gradient 3')

        call self%add_model_field('Reynolds-13 - Gradient 1')
        call self%add_model_field('Reynolds-13 - Gradient 2')
        call self%add_model_field('Reynolds-13 - Gradient 3')

        call self%add_model_field('Reynolds-23 - Gradient 1')
        call self%add_model_field('Reynolds-23 - Gradient 2')
        call self%add_model_field('Reynolds-23 - Gradient 3')

        call self%add_model_field('Anisotropy-11')
        call self%add_model_field('Anisotropy-22')
        call self%add_model_field('Anisotropy-33')
        call self%add_model_field('Anisotropy-12')
        call self%add_model_field('Anisotropy-13')
        call self%add_model_field('Anisotropy-23')


        ! k
        call self%add_model_field('Turbulence Kinetic Energy')
 
        ! gradient k
        call self%add_model_field('Turbulence Kinetic Energy - Gradient 1')
        call self%add_model_field('Turbulence Kinetic Energy - Gradient 2')
        call self%add_model_field('Turbulence Kinetic Energy - Gradient 3')


        ! epsilon
        call self%add_model_field('Turbulence Isotropic Dissipation Rate')

        ! (density/omega)*max(dk/dx_{k}*domega/dx_{k}, 0)
        call self%add_model_field('Omega Source Term')

        call self%add_model_field('Equivalent Eddy Viscosity')
        !call self%add_model_field('Turbulent Viscosity')
        !call self%add_model_field('Second Coefficient of Turbulent Viscosity')
        !call self%add_model_field('Turbulent Thermal Conductivity')
!
        ! Check if input from 'models.nml' is available.
        !   1: if available, read and set self%mu
        !   2: if not available, do nothing and mu retains default value
        !
        inquire(file='models.nml', exist=file_exists)
        if (file_exists) then
            open(newunit=unit,form='formatted',file='models.nml')
            read(unit,nml=blended_rstm,iostat=msg)
            if (msg == 0) self%blend = blend 
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
        class(rstm_ssglrrw_blended_turbulence_quantities_t),     intent(in)      :: self
        type(chidg_worker_t),   intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable ::              &
            density, mu,    &
            grad1_density,  grad2_density,  grad3_density,              &
            density_omega, omega,                                       &
            grad1_density_omega, grad2_density_omega, grad3_density_omega, &
            grad1_omega,    grad2_omega,    grad3_omega,                    &
            reynolds_11, reynolds_22, reynolds_33,                      &
            reynolds_12, reynolds_13, reynolds_23,                      &
            grad1_reynolds_11, grad2_reynolds_11, grad3_reynolds_11,    &
            grad1_reynolds_22, grad2_reynolds_22, grad3_reynolds_22,    &
            grad1_reynolds_33, grad2_reynolds_33, grad3_reynolds_33,    &
            grad1_reynolds_12, grad2_reynolds_12, grad3_reynolds_12,    &
            grad1_reynolds_13, grad2_reynolds_13, grad3_reynolds_13,    &
            grad1_reynolds_23, grad2_reynolds_23, grad3_reynolds_23,    &
            grad1_density_reynolds_11, grad2_density_reynolds_11, grad3_density_reynolds_11,    &
            grad1_density_reynolds_22, grad2_density_reynolds_22, grad3_density_reynolds_22,    &
            grad1_density_reynolds_33, grad2_density_reynolds_33, grad3_density_reynolds_33,    &
            grad1_density_reynolds_12, grad2_density_reynolds_12, grad3_density_reynolds_12,    &
            grad1_density_reynolds_13, grad2_density_reynolds_13, grad3_density_reynolds_13,    &
            grad1_density_reynolds_23, grad2_density_reynolds_23, grad3_density_reynolds_23,    &
            anisotropy_11, anisotropy_22, anisotropy_33,                      &
            anisotropy_12, anisotropy_13, anisotropy_23,                      &
            grad1_k_t, grad2_k_t, grad3_k_t,                            &
            mu_t, k_t, epsilon_t, temp1, temp2, omega_source_term, invdensity, grad_omega_sq

        real(rk), allocatable :: det_r(:)
        !
        ! Interpolate solution to quadrature nodes
        !
        density = worker%get_field('Density',    'value')
        grad1_density = worker%get_field('Density', 'grad1')
        grad2_density = worker%get_field('Density', 'grad2')
        grad3_density = worker%get_field('Density', 'grad3')

        
        if (any(ieee_is_nan(density(:)%x_ad_))) print *, 'turb quant density is nan'
        
        invdensity = ONE/density
        density_omega = worker%get_field('Density * Omega',    'value')
        grad1_density_omega = worker%get_field('Density * Omega', 'grad1')
        grad2_density_omega = worker%get_field('Density * Omega', 'grad2')
        grad3_density_omega = worker%get_field('Density * Omega', 'grad3')
        omega = (density_omega*invdensity)
        
        !if (any(omega%x_ad_<0.0_rk)) print *, 'warning omega < 0, ', worker%interpolation_source
        grad1_omega = invdensity*(grad1_density_omega-omega*grad1_density)
        grad2_omega = invdensity*(grad2_density_omega-omega*grad2_density)
        grad3_omega = invdensity*(grad3_density_omega-omega*grad3_density)

        grad_omega_sq = grad1_omega*grad1_omega + grad2_omega*grad2_omega + grad3_omega*grad3_omega

        call worker%store_model_field('Omega', 'value', (omega))
        call worker%store_model_field('Omega - Gradient 1', 'value', grad1_omega)
        call worker%store_model_field('Omega - Gradient 2', 'value', grad2_omega)
        call worker%store_model_field('Omega - Gradient 3', 'value', grad3_omega)
        call worker%store_model_field('Omega Gradient Squared', 'value', grad_omega_sq)

        reynolds_11 = worker%get_field('Density * Reynolds-11', 'value')
        reynolds_22 = worker%get_field('Density * Reynolds-22', 'value')
        reynolds_33 = worker%get_field('Density * Reynolds-33', 'value')
        reynolds_12 = worker%get_field('Density * Reynolds-12', 'value')
        reynolds_13 = worker%get_field('Density * Reynolds-13', 'value')
        reynolds_23 = worker%get_field('Density * Reynolds-23', 'value')

        !if ((any(reynolds_11(:)%x_ad_ < 0.0_rk)) .or. (any(reynolds_22(:)%x_ad_ < 0.0_rk)) .or. (any(reynolds_33(:)%x_ad_ < 0.0_rk))) then
        !    print *, 'warning, negative R diag, ', worker%interpolation_source
        !end if
        !
        !if (any(abs(reynolds_12(:)%x_ad_) > sqrt(abs(reynolds_11(:)%x_ad_*reynolds_22(:)%x_ad_)))) print *, 'warning, r_12 unrealizable, ', worker%interpolation_source
        !if (any(abs(reynolds_13(:)%x_ad_) > sqrt(abs(reynolds_11(:)%x_ad_*reynolds_33(:)%x_ad_)))) print *, 'warning, r_13 unrealizable, ', worker%interpolation_source
        !if (any(abs(reynolds_23(:)%x_ad_) > sqrt(abs(reynolds_22(:)%x_ad_*reynolds_33(:)%x_ad_)))) print *, 'warning, r_23 unrealizable, ', worker%interpolation_source

        !det_r =  reynolds_11(:)%x_ad_*(reynolds_22(:)%x_ad_*reynolds_33(:)%x_ad_-reynolds_23(:)%x_ad_*reynolds_23(:)%x_ad_) &
        !        -reynolds_12(:)%x_ad_*(reynolds_12(:)%x_ad_*reynolds_33(:)%x_ad_-reynolds_23(:)%x_ad_*reynolds_13(:)%x_ad_) &
        !        +reynolds_13(:)%x_ad_*(reynolds_12(:)%x_ad_*reynolds_23(:)%x_ad_-reynolds_22(:)%x_ad_*reynolds_13(:)%x_ad_) 

        !if (any(det_r<0.0_rk)) print *, 'warning det R < 0, ', worker%interpolation_source

        !!call worker%store_model_field('Reynolds-Stress-11', 'value', -(reynolds_11))
        !!call worker%store_model_field('Reynolds-Stress-22', 'value', -(reynolds_22))
        !!call worker%store_model_field('Reynolds-Stress-33', 'value', -(reynolds_33))
        !!call worker%store_model_field('Reynolds-Stress-12', 'value', -reynolds_12)
        !!call worker%store_model_field('Reynolds-Stress-13', 'value', -reynolds_13)
        !!call worker%store_model_field('Reynolds-Stress-23', 'value', -reynolds_23)


        reynolds_11 = invdensity*reynolds_11
        reynolds_22 = invdensity*reynolds_22
        reynolds_33 = invdensity*reynolds_33
        reynolds_12 = invdensity*reynolds_12
        reynolds_13 = invdensity*reynolds_13
        reynolds_23 = invdensity*reynolds_23


        !call worker%store_model_field('Reynolds-11', 'value', (reynolds_11))
        !call worker%store_model_field('Reynolds-22', 'value', (reynolds_22))
        !call worker%store_model_field('Reynolds-33', 'value', (reynolds_33))
        !call worker%store_model_field('Reynolds-12', 'value', reynolds_12)
        !call worker%store_model_field('Reynolds-13', 'value', reynolds_13)
        !call worker%store_model_field('Reynolds-23', 'value', reynolds_23)

        grad1_density_reynolds_11 = worker%get_field('Density * Reynolds-11', 'grad1')
        grad2_density_reynolds_11 = worker%get_field('Density * Reynolds-11', 'grad2')
        grad3_density_reynolds_11 = worker%get_field('Density * Reynolds-11', 'grad3')

        grad1_density_reynolds_22 = worker%get_field('Density * Reynolds-22', 'grad1')
        grad2_density_reynolds_22 = worker%get_field('Density * Reynolds-22', 'grad2')
        grad3_density_reynolds_22 = worker%get_field('Density * Reynolds-22', 'grad3')

        grad1_density_reynolds_33 = worker%get_field('Density * Reynolds-33', 'grad1')
        grad2_density_reynolds_33 = worker%get_field('Density * Reynolds-33', 'grad2')
        grad3_density_reynolds_33 = worker%get_field('Density * Reynolds-33', 'grad3')

        grad1_density_reynolds_12 = worker%get_field('Density * Reynolds-12', 'grad1')
        grad2_density_reynolds_12 = worker%get_field('Density * Reynolds-12', 'grad2')
        grad3_density_reynolds_12 = worker%get_field('Density * Reynolds-12', 'grad3')

        grad1_density_reynolds_13 = worker%get_field('Density * Reynolds-13', 'grad1')
        grad2_density_reynolds_13 = worker%get_field('Density * Reynolds-13', 'grad2')
        grad3_density_reynolds_13 = worker%get_field('Density * Reynolds-13', 'grad3')

        grad1_density_reynolds_23 = worker%get_field('Density * Reynolds-23', 'grad1')
        grad2_density_reynolds_23 = worker%get_field('Density * Reynolds-23', 'grad2')
        grad3_density_reynolds_23 = worker%get_field('Density * Reynolds-23', 'grad3')


        grad1_reynolds_11 = invdensity*(grad1_density_reynolds_11 - reynolds_11*grad1_density)
        grad2_reynolds_11 = invdensity*(grad2_density_reynolds_11 - reynolds_11*grad2_density)
        grad3_reynolds_11 = invdensity*(grad3_density_reynolds_11 - reynolds_11*grad3_density)

        grad1_reynolds_22 = invdensity*(grad1_density_reynolds_22 - reynolds_22*grad1_density)
        grad2_reynolds_22 = invdensity*(grad2_density_reynolds_22 - reynolds_22*grad2_density)
        grad3_reynolds_22 = invdensity*(grad3_density_reynolds_22 - reynolds_22*grad3_density)

        grad1_reynolds_33 = invdensity*(grad1_density_reynolds_33 - reynolds_33*grad1_density)
        grad2_reynolds_33 = invdensity*(grad2_density_reynolds_33 - reynolds_33*grad2_density)
        grad3_reynolds_33 = invdensity*(grad3_density_reynolds_33 - reynolds_33*grad3_density)
        
        grad1_reynolds_12 = invdensity*(grad1_density_reynolds_12 - reynolds_12*grad1_density)
        grad2_reynolds_12 = invdensity*(grad2_density_reynolds_12 - reynolds_12*grad2_density)
        grad3_reynolds_12 = invdensity*(grad3_density_reynolds_12 - reynolds_12*grad3_density)

        grad1_reynolds_13 = invdensity*(grad1_density_reynolds_13 - reynolds_13*grad1_density)
        grad2_reynolds_13 = invdensity*(grad2_density_reynolds_13 - reynolds_13*grad2_density)
        grad3_reynolds_13 = invdensity*(grad3_density_reynolds_13 - reynolds_13*grad3_density)

        grad1_reynolds_23 = invdensity*(grad1_density_reynolds_23 - reynolds_23*grad1_density)
        grad2_reynolds_23 = invdensity*(grad2_density_reynolds_23 - reynolds_23*grad2_density)
        grad3_reynolds_23 = invdensity*(grad3_density_reynolds_23 - reynolds_23*grad3_density)


        call worker%store_model_field('Reynolds-11 - Gradient 1', 'value', grad1_reynolds_11)
        call worker%store_model_field('Reynolds-11 - Gradient 2', 'value', grad2_reynolds_11)
        call worker%store_model_field('Reynolds-11 - Gradient 3', 'value', grad3_reynolds_11)

        call worker%store_model_field('Reynolds-22 - Gradient 1', 'value', grad1_reynolds_22)
        call worker%store_model_field('Reynolds-22 - Gradient 2', 'value', grad2_reynolds_22)
        call worker%store_model_field('Reynolds-22 - Gradient 3', 'value', grad3_reynolds_22)

        call worker%store_model_field('Reynolds-33 - Gradient 1', 'value', grad1_reynolds_33)
        call worker%store_model_field('Reynolds-33 - Gradient 2', 'value', grad2_reynolds_33)
        call worker%store_model_field('Reynolds-33 - Gradient 3', 'value', grad3_reynolds_33)

        call worker%store_model_field('Reynolds-12 - Gradient 1', 'value', grad1_reynolds_12)
        call worker%store_model_field('Reynolds-12 - Gradient 2', 'value', grad2_reynolds_12)
        call worker%store_model_field('Reynolds-12 - Gradient 3', 'value', grad3_reynolds_12)

        call worker%store_model_field('Reynolds-13 - Gradient 1', 'value', grad1_reynolds_13)
        call worker%store_model_field('Reynolds-13 - Gradient 2', 'value', grad2_reynolds_13)
        call worker%store_model_field('Reynolds-13 - Gradient 3', 'value', grad3_reynolds_13)

        call worker%store_model_field('Reynolds-23 - Gradient 1', 'value', grad1_reynolds_23)
        call worker%store_model_field('Reynolds-23 - Gradient 2', 'value', grad2_reynolds_23)
        call worker%store_model_field('Reynolds-23 - Gradient 3', 'value', grad3_reynolds_23)

        !
        ! Get realizable Reynolds stress tensor
        !
        reynolds_11 = worker%get_field('Reynolds-11', 'value')
        reynolds_22 = worker%get_field('Reynolds-22', 'value')
        reynolds_33 = worker%get_field('Reynolds-33', 'value')
        reynolds_12 = worker%get_field('Reynolds-12', 'value')
        reynolds_13 = worker%get_field('Reynolds-13', 'value')
        reynolds_23 = worker%get_field('Reynolds-23', 'value')


        k_t = 0.5_rk*((reynolds_11)+(reynolds_22)+(reynolds_33))
        grad1_k_t = 0.5_rk*(grad1_reynolds_11 + grad1_reynolds_22 + grad1_reynolds_33)
        grad2_k_t = 0.5_rk*(grad2_reynolds_11 + grad2_reynolds_22 + grad2_reynolds_33)
        grad3_k_t = 0.5_rk*(grad3_reynolds_11 + grad3_reynolds_22 + grad3_reynolds_33)

        call worker%store_model_field('Turbulence Kinetic Energy', 'value', (k_t))
        call worker%store_model_field('Turbulence Kinetic Energy - Gradient 1', 'value', grad1_k_t)
        call worker%store_model_field('Turbulence Kinetic Energy - Gradient 2', 'value', grad2_k_t)
        call worker%store_model_field('Turbulence Kinetic Energy - Gradient 3', 'value', grad3_k_t)


        epsilon_t = SSG_LRRW_cmu*k_t*exp(omega)
        call worker%store_model_field('Turbulence Isotropic Dissipation Rate', 'value', (epsilon_t))

        anisotropy_11 = reynolds_11/(k_t+1.0e-11_rk)-(2.0_rk/3.0_rk)
        anisotropy_22 = reynolds_22/(k_t+1.0e-11_rk)-(2.0_rk/3.0_rk)
        anisotropy_33 = reynolds_33/(k_t+1.0e-11_rk)-(2.0_rk/3.0_rk)
        anisotropy_12 = reynolds_12/(k_t+1.0e-11_rk)
        anisotropy_13 = reynolds_13/(k_t+1.0e-11_rk)
        anisotropy_23 = reynolds_23/(k_t+1.0e-11_rk)
        
        call worker%store_model_field('Anisotropy-11', 'value', anisotropy_11)
        call worker%store_model_field('Anisotropy-22', 'value', anisotropy_22)
        call worker%store_model_field('Anisotropy-33', 'value', anisotropy_33)
        call worker%store_model_field('Anisotropy-12', 'value', anisotropy_12)
        call worker%store_model_field('Anisotropy-13', 'value', anisotropy_13)
        call worker%store_model_field('Anisotropy-23', 'value', anisotropy_23)

        temp1 = grad1_k_t*grad1_omega+grad2_k_t*grad2_omega+grad3_k_t*grad3_omega
        temp2 = temp1*sin_ramp(temp1, 0.0_rk, 10.0_rk*rstm_ssglrrw_k_infty)
        !temp2 = 0.5_rk*(abs(temp1) + temp1)
        !temp2 = max(temp1, ZERO)
        omega_source_term = (density*exp(-omega))*temp2

        call worker%store_model_field('Omega Source Term', 'value', omega_source_term)


        mu_t = density*k_t*exp(-omega)
        call worker%store_model_field('Equivalent Eddy Viscosity', 'value', (mu_t))

        !call worker%store_model_field('Turbulent Viscosity', 'value', ZERO*(mu_t))
        !call worker%store_model_field('Second Coefficient of Turbulent Viscosity', 'value', ZERO*(mu_t))
        !call worker%store_model_field('Turbulent Thermal Conductivity', 'value', cp*(mu_t)/0.9_rk)
    end subroutine compute
    !***************************************************************************************




end module model_rstm_ssglrrw_blended_turbulence_quantities
