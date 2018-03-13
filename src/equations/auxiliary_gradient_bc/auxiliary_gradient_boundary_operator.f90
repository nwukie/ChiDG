module auxiliary_gradient_boundary_operator
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO,ONE,TWO,HALF
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    use ieee_arithmetic
    implicit none

    private



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/5/2018
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: auxiliary_gradient_boundary_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type auxiliary_gradient_boundary_operator_t
    !********************************************************************************



contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/5/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(auxiliary_gradient_boundary_operator_t),   intent(inout) :: self
        
        ! Set operator name
        call self%set_name("Auxiliary Gradient Boundary Average Operator")

        ! Set operator type
        call self%set_operator_type("Boundary Diffusive Operator")

        ! Set operator equations
        call self%add_primary_field('Density_TEMP')
        call self%add_primary_field('Velocity-1_TEMP')
        call self%add_primary_field('Pressure_TEMP')

    end subroutine init
    !********************************************************************************







    !>  Compute the diffusive boundary flux for scalar linear diffusion.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      mesh    Mesh data
    !!  @param[inout]   sdata   Solver data. Solution, RHS, Linearization etc.
    !!  @param[in]      ielem   Element index
    !!  @param[in]      iface   Face index
    !!  @param[in]      iblk    Block index indicating the linearization direction
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(auxiliary_gradient_boundary_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:)   ::                      &
            grad1_p_bc_m,       grad2_p_bc_m,       grad3_p_bc_m,       &
            grad1_p_bc_p,       grad2_p_bc_p,       grad3_p_bc_p,       &
            grad1_p_m,          grad2_p_m,          grad3_p_m,          &
            grad1_p_p,          grad2_p_p,          grad3_p_p,          &
            grad1_density_bc_m, grad2_density_bc_m, grad3_density_bc_m, &
            grad1_density_bc_p, grad2_density_bc_p, grad3_density_bc_p, &
            grad1_density_m,    grad2_density_m,    grad3_density_m,    &
            grad1_density_p,    grad2_density_p,    grad3_density_p,    &
            grad1_vel1_bc_m,    grad2_vel1_bc_m,    grad3_vel1_bc_m,    &
            grad1_vel1_bc_p,    grad2_vel1_bc_p,    grad3_vel1_bc_p,    &
            grad1_vel1_m,       grad2_vel1_m,       grad3_vel1_m,       &
            grad1_vel1_p,       grad2_vel1_p,       grad3_vel1_p,       &
            flux_1_m, flux_2_m, flux_3_m,                               &
            flux_1_p, flux_2_p, flux_3_p


        

        !-----------------------------------
        !   Auxiliary Density
        !-----------------------------------
        ! density_bc
        grad1_density_bc_m = worker%get_field('Density_TEMP', 'grad1', 'face interior')
        grad2_density_bc_m = worker%get_field('Density_TEMP', 'grad2', 'face interior')
        grad3_density_bc_m = worker%get_field('Density_TEMP', 'grad3', 'face interior')

        grad1_density_bc_p = worker%get_field('Density_TEMP', 'grad1', 'face exterior')
        grad2_density_bc_p = worker%get_field('Density_TEMP', 'grad2', 'face exterior')
        grad3_density_bc_p = worker%get_field('Density_TEMP', 'grad3', 'face exterior')

        ! SIGMA
        grad1_density_m   = worker%get_field('Density', 'grad1', 'face interior', override_lift=.true.)
        grad2_density_m   = worker%get_field('Density', 'grad2', 'face interior', override_lift=.true.)
        grad3_density_m   = worker%get_field('Density', 'grad3', 'face interior', override_lift=.true.)

        grad1_density_p   = worker%get_field('Density', 'grad1', 'face exterior', override_lift=.true.)
        grad2_density_p   = worker%get_field('Density', 'grad2', 'face exterior', override_lift=.true.)
        grad3_density_p   = worker%get_field('Density', 'grad3', 'face exterior', override_lift=.true.)


        ! Compute flux from each side
        flux_1_m = grad1_density_bc_m - grad1_density_m
        flux_2_m = grad2_density_bc_m - grad2_density_m
        flux_3_m = grad3_density_bc_m - grad3_density_m

        flux_1_p = grad1_density_bc_p - grad1_density_p
        flux_2_p = grad2_density_bc_p - grad2_density_p
        flux_3_p = grad3_density_bc_p - grad3_density_p


        ! Integrate flux
        call worker%integrate_boundary_average('Density_TEMP','Diffusion',      &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)




        !-----------------------------------
        !   Auxiliary Velocity-1
        !-----------------------------------
        ! vel1_bc
        grad1_vel1_bc_m = worker%get_field('Velocity-1_TEMP', 'grad1', 'face interior')
        grad2_vel1_bc_m = worker%get_field('Velocity-1_TEMP', 'grad2', 'face interior')
        grad3_vel1_bc_m = worker%get_field('Velocity-1_TEMP', 'grad3', 'face interior')

        grad1_vel1_bc_p = worker%get_field('Velocity-1_TEMP', 'grad1', 'face exterior')
        grad2_vel1_bc_p = worker%get_field('Velocity-1_TEMP', 'grad2', 'face exterior')
        grad3_vel1_bc_p = worker%get_field('Velocity-1_TEMP', 'grad3', 'face exterior')

        ! SIGMA
        grad1_vel1_m   = worker%get_field('Velocity-1 Gradient-1', 'value', 'face interior')
        grad2_vel1_m   = worker%get_field('Velocity-1 Gradient-2', 'value', 'face interior')
        grad3_vel1_m   = worker%get_field('Velocity-1 Gradient-3', 'value', 'face interior')

        grad1_vel1_p   = worker%get_field('Velocity-1 Gradient-1', 'value', 'face exterior')
        grad2_vel1_p   = worker%get_field('Velocity-1 Gradient-2', 'value', 'face exterior')
        grad3_vel1_p   = worker%get_field('Velocity-1 Gradient-3', 'value', 'face exterior')


        ! Compute flux from each side
        flux_1_m = grad1_vel1_bc_m - grad1_vel1_m
        flux_2_m = grad2_vel1_bc_m - grad2_vel1_m
        flux_3_m = grad3_vel1_bc_m - grad3_vel1_m

        flux_1_p = grad1_vel1_bc_p - grad1_vel1_p
        flux_2_p = grad2_vel1_bc_p - grad2_vel1_p
        flux_3_p = grad3_vel1_bc_p - grad3_vel1_p


        ! Integrate flux
        call worker%integrate_boundary_average('Velocity-1_TEMP','Diffusion',   &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)



        !-----------------------------------
        !   Auxiliary Pressure
        !-----------------------------------
        ! P_bc
        grad1_p_bc_m = worker%get_field('Pressure_TEMP', 'grad1', 'face interior')
        grad2_p_bc_m = worker%get_field('Pressure_TEMP', 'grad2', 'face interior')
        grad3_p_bc_m = worker%get_field('Pressure_TEMP', 'grad3', 'face interior')

        grad1_p_bc_p = worker%get_field('Pressure_TEMP', 'grad1', 'face exterior')
        grad2_p_bc_p = worker%get_field('Pressure_TEMP', 'grad2', 'face exterior')
        grad3_p_bc_p = worker%get_field('Pressure_TEMP', 'grad3', 'face exterior')

        ! SIGMA
        grad1_p_m   = worker%get_field('Pressure Gradient - 1', 'value', 'face interior')
        grad2_p_m   = worker%get_field('Pressure Gradient - 2', 'value', 'face interior')
        grad3_p_m   = worker%get_field('Pressure Gradient - 3', 'value', 'face interior')

        grad1_p_p   = worker%get_field('Pressure Gradient - 1', 'value', 'face exterior')
        grad2_p_p   = worker%get_field('Pressure Gradient - 2', 'value', 'face exterior')
        grad3_p_p   = worker%get_field('Pressure Gradient - 3', 'value', 'face exterior')


        ! Compute flux from each side
        flux_1_m = grad1_p_bc_m - grad1_p_m
        flux_2_m = grad2_p_bc_m - grad2_p_m
        flux_3_m = grad3_p_bc_m - grad3_p_m

        flux_1_p = grad1_p_bc_p - grad1_p_p
        flux_2_p = grad2_p_bc_p - grad2_p_p
        flux_3_p = grad3_p_bc_p - grad3_p_p


        ! Integrate flux
        call worker%integrate_boundary_average('Pressure_TEMP','Diffusion',     &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)



    end subroutine compute
    !**************************************************************************************************




end module auxiliary_gradient_boundary_operator
