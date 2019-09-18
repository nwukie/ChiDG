module model_pde_smoothed_artificial_viscosity
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: THREE, TWO, ONE, ZERO, HALF, PI, THIRD
    use mod_fluid
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use DNAD_D
    use ieee_arithmetic
    implicit none


    


    !>  Artificial viscosity model based on
    !!  A physics-based shock capturing method for unsteady laminar and turbulent flows
    !!      Fernandez et al, 2018, AIAA SciTech Forum
    !!
    !!  Model Fields:
    !!      - Artifical Bulk Viscosity
    !!      - Artifical Shear Viscosity
    !!      - Artifical Thermal Conductivity 
    !!
    !!  @author Eric M. Wolf
    !!  @date   07/11/2018
    !!
    !---------------------------------------------------------------------------------------
    type, extends(model_t)  :: pde_smoothed_artificial_viscosity_t


    contains

        procedure   :: init
        procedure   :: compute

    end type pde_smoothed_artificial_viscosity_t
    !***************************************************************************************





contains




    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    07/11/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)   
        class(pde_smoothed_artificial_viscosity_t), intent(inout)   :: self
        




        call self%set_name('PDE Smoothed Artificial Viscosity')
        call self%set_dependency('f(Grad(Q))')

        call self%add_model_field('Smoothed Artificial Bulk Viscosity')
        call self%add_model_field('Smoothed Artificial Shear Viscosity')
        call self%add_model_field('Smoothed Artificial Thermal Conductivity')

        

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
        class(pde_smoothed_artificial_viscosity_t),   intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            density, vel1, vel2, vel3, T, c, wave_speed, shock_sensor, &
            art_bulk_viscosity, art_shear_viscosity, art_thermal_conductivity, eps, lambda

        type(AD_D) :: theta_l, theta_h
        
        integer(ik)             :: order, inode
        real(rk)                :: h(3), hbar, lambda_max, theta_l_r, theta_h_r
        real(rk)                :: Pr_star  = 0.9_rk     

        eps =  worker%get_field('Artificial Viscosity', 'value')
        eps = 5.0_rk*eps
        art_bulk_viscosity          =  worker%get_field('Artificial Viscosity', 'value')

        lambda = worker%get_field('Maximum Wave Speed', 'value')

        lambda_max = maxval(lambda(:)%x_ad_)


        order = worker%solution_order('interior')
        if (order==0) order = 1
        h = worker%element_size('interior')

        hbar = THIRD*(h(1)+h(2)+h(3))
        !hbar = h(1)

        do inode = 1, size(art_bulk_viscosity)
            theta_l = 0.01_rk*5.0_rk*lambda(inode)*hbar/real(order,rk)
            theta_h = 5.0_rk*lambda(inode)*hbar/real(order,rk)

            theta_l_r = 0.01_rk*lambda_max*hbar/real(order,rk)
            theta_h_r = lambda_max*hbar/real(order,rk)

            if (eps(inode) < theta_l) then
                art_bulk_viscosity(inode) = ZERO
            else if ((theta_l <= eps(inode)) .and. (eps(inode) < theta_h)) then
                art_bulk_viscosity(inode) = HALF*theta_h*(sin(PI*((eps(inode)-theta_l)/(theta_h-theta_l)-HALF))+ONE)
            else 
                art_bulk_viscosity(inode) = theta_h
            end if

!            if (eps(inode) < theta_l_r) then
!                art_bulk_viscosity(inode) = ZERO
!            else if ((theta_l_r <= eps(inode)) .and. (eps(inode) < theta_h_r)) then
!                art_bulk_viscosity(inode) = HALF*theta_h_r*(sin(PI*((eps(inode)-theta_l_r)/(theta_h_r-theta_l_r)-HALF))+ONE)
!            else 
!                art_bulk_viscosity(inode) = theta_h_r
!            end if

            art_bulk_viscosity(inode)%xp_ad_(:) = ZERO

        end do

        art_shear_viscosity         = art_bulk_viscosity
        art_thermal_conductivity = (cp/Pr_star*art_bulk_viscosity + cp*art_bulk_viscosity)




        !
        ! Contribute laminar viscosity
        !
        call worker%store_model_field('Smoothed Artificial Bulk Viscosity', 'value', art_bulk_viscosity)
        call worker%store_model_field('Smoothed Artificial Shear Viscosity', 'value', art_shear_viscosity)
        call worker%store_model_field('Smoothed Artificial Thermal Conductivity', 'value', art_thermal_conductivity)


    end subroutine compute
    !***************************************************************************************




end module model_pde_smoothed_artificial_viscosity
