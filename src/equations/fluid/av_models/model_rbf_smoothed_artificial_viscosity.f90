module model_rbf_smoothed_artificial_viscosity
#include <messenger.h>
    use mod_kinds,          only: rk
    use mod_constants,      only: THREE, TWO, ONE, ZERO
    use mod_fluid
    use type_model,         only: model_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_radial_basis_function
    use mod_radial_basis_function
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
    type, extends(model_t)  :: rbf_smoothed_artificial_viscosity_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rbf_smoothed_artificial_viscosity_t
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
        class(rbf_smoothed_artificial_viscosity_t), intent(inout)   :: self

        call self%set_name('RBF Smoothed Artificial Viscosity')
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
        class(rbf_smoothed_artificial_viscosity_t),   intent(in)      :: self
        type(chidg_worker_t),       intent(inout)   :: worker

        type(AD_D), dimension(:),   allocatable :: &
            density, vel1, vel2, vel3, T, c, wave_speed, shock_sensor, &
            art_bulk_viscosity, art_shear_viscosity, art_thermal_conductivity

        integer(ik)             :: p, ii, nrbf, inode
        real(rk), allocatable   :: eval_node1(:), eval_node2(:), eval_node3(:)
        real(rk)                :: eval_node(3), center(3), radius(3)
        class(radial_basis_function_t), allocatable  :: rbf
        
        call create_radial_basis_function(rbf,'wc2_rbf')

        eval_node1 = worker%coordinate('1')
        eval_node2 = worker%coordinate('2')
        eval_node3 = worker%coordinate('3')
        
        density = worker%get_field('Density','value')

        art_bulk_viscosity = ZERO*density
        art_shear_viscosity = ZERO*density
        art_thermal_conductivity = ZERO*density
        nrbf = size(worker%solverdata%rbf_center(:,1))
        do ii = 1, nrbf
            center = worker%solverdata%rbf_center(ii,:)
            radius = 2.5_rk*sqrt(3.0_rk)*worker%solverdata%rbf_radius(ii,:)

            radius(2) = 1.0e12_rk
            radius(3) = 1.0e12_rk

            !print *, 'rbf center', center
            !print *, 'rbf radius', radius
            !print *, 'rbf asv', worker%solverdata%rbf_artificial_shear_viscosity(ii)
            
            if (ieee_is_nan(worker%solverdata%rbf_artificial_bulk_viscosity(ii))) print *, 'art bulk visc is nan'
            if (ieee_is_nan(worker%solverdata%rbf_artificial_shear_viscosity(ii))) print *, 'art shear visc is nan'
            if (ieee_is_nan(worker%solverdata%rbf_artificial_thermal_conductivity(ii))) print *, 'art therm cond is nan'
           
            do inode = 1, size(eval_node1)
                eval_node(1) = eval_node1(inode)
                eval_node(2) = eval_node2(inode)
                eval_node(3) = eval_node3(inode)
            
                art_bulk_viscosity(inode) = art_bulk_viscosity(inode) + worker%solverdata%rbf_artificial_bulk_viscosity(ii)*rbf%compute(eval_node, center, radius)
                art_shear_viscosity(inode) = art_shear_viscosity(inode) + worker%solverdata%rbf_artificial_shear_viscosity(ii)*rbf%compute(eval_node, center, radius)
                art_thermal_conductivity(inode) = art_thermal_conductivity(inode) + worker%solverdata%rbf_artificial_thermal_conductivity(ii)*rbf%compute(eval_node, center, radius)
            end do
            
        end do

       
        !
        ! Contribute laminar viscosity
        !
        call worker%store_model_field('Smoothed Artificial Bulk Viscosity', 'value', art_bulk_viscosity)
        call worker%store_model_field('Smoothed Artificial Shear Viscosity', 'value', art_shear_viscosity)
        call worker%store_model_field('Smoothed Artificial Thermal Conductivity', 'value', art_thermal_conductivity)


    end subroutine compute
    !***************************************************************************************




end module model_rbf_smoothed_artificial_viscosity
