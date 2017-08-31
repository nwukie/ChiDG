module MMD_boundary_operator
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: ZERO,ONE,TWO,HALF
    use type_operator,              only: operator_t
    use type_chidg_worker,          only: chidg_worker_t
    use type_properties,            only: properties_t
    use DNAD_D
    use ieee_arithmetic
    implicit none

    private



    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: MMD_boundary_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type MMD_boundary_operator_t
    !********************************************************************************



contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(MMD_boundary_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Mesh Motion Diffusion Boundary Average Operator")

        !
        ! Set operator type
        !
        call self%set_operator_type("Boundary Diffusive Operator")

        !
        ! Set operator equations
        !
        call self%add_primary_field("grid_displacement1")
        call self%add_primary_field("grid_displacement2")
        call self%add_primary_field("grid_displacement3")

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
        class(MMD_boundary_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop



        type(AD_D), allocatable, dimension(:)   ::  &
            grad1_u1_m, grad2_u1_m, grad3_u1_m,     &
            grad1_u1_p, grad2_u1_p, grad3_u1_p,     &
            grad1_u2_m, grad2_u2_m, grad3_u2_m,     &
            grad1_u2_p, grad2_u2_p, grad3_u2_p,     &
            grad1_u3_m, grad2_u3_m, grad3_u3_m,     &
            grad1_u3_p, grad2_u3_p, grad3_u3_p,     &
            flux_1_m, flux_2_m, flux_3_m,           &
            flux_1_p, flux_2_p, flux_3_p,           &
            mu_m, mu_p



        !
        ! Interpolate solution to quadrature nodes
        !
        grad1_u1_m = worker%get_field('grid_displacement1', 'grad1', 'face interior')
        grad2_u1_m = worker%get_field('grid_displacement1', 'grad2', 'face interior')
        grad3_u1_m = worker%get_field('grid_displacement1', 'grad3', 'face interior')


        grad1_u1_p = worker%get_field('grid_displacement1', 'grad1', 'face exterior')
        grad2_u1_p = worker%get_field('grid_displacement1', 'grad2', 'face exterior')
        grad3_u1_p = worker%get_field('grid_displacement1', 'grad3', 'face exterior')

        grad1_u2_m = worker%get_field('grid_displacement2', 'grad1', 'face interior')
        grad2_u2_m = worker%get_field('grid_displacement2', 'grad2', 'face interior')
        grad3_u2_m = worker%get_field('grid_displacement2', 'grad3', 'face interior')


        grad1_u2_p = worker%get_field('grid_displacement2', 'grad1', 'face exterior')
        grad2_u2_p = worker%get_field('grid_displacement2', 'grad2', 'face exterior')
        grad3_u2_p = worker%get_field('grid_displacement2', 'grad3', 'face exterior')

        grad1_u3_m = worker%get_field('grid_displacement3', 'grad1', 'face interior')
        grad2_u3_m = worker%get_field('grid_displacement3', 'grad2', 'face interior')
        grad3_u3_m = worker%get_field('grid_displacement3', 'grad3', 'face interior')


        grad1_u3_p = worker%get_field('grid_displacement3', 'grad1', 'face exterior')
        grad2_u3_p = worker%get_field('grid_displacement3', 'grad2', 'face exterior')
        grad3_u3_p = worker%get_field('grid_displacement3', 'grad3', 'face exterior')






        !
        ! Compute scalar coefficient
        !
        mu_m = worker%get_field('Mesh Motion Diffusion Coefficient', 'value', 'face interior')
        mu_p = worker%get_field('Mesh Motion Diffusion Coefficient', 'value', 'face exterior')




        !======================================
        !               GD1
        !======================================
        flux_1_m = -mu_m*grad1_u1_m
        flux_2_m = -mu_m*grad2_u1_m
        flux_3_m = -mu_m*grad3_u1_m

        flux_1_p = -mu_p*grad1_u1_p
        flux_2_p = -mu_p*grad2_u1_p
        flux_3_p = -mu_p*grad3_u1_p


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_average('grid_displacement1','Diffusion',    &
                                                flux_1_m, flux_2_m, flux_3_m,       &
                                                flux_1_p, flux_2_p, flux_3_p)


        !======================================
        !               GD2
        !======================================
        flux_1_m = -mu_m*grad1_u2_m
        flux_2_m = -mu_m*grad2_u2_m
        flux_3_m = -mu_m*grad3_u2_m

        flux_1_p = -mu_p*grad1_u2_p
        flux_2_p = -mu_p*grad2_u2_p
        flux_3_p = -mu_p*grad3_u2_p

        !
        ! Integrate flux
        !
        call worker%integrate_boundary_average('grid_displacement2','Diffusion',    &
                                                flux_1_m, flux_2_m, flux_3_m,       &
                                                flux_1_p, flux_2_p, flux_3_p)


        !======================================
        !               GD3
        !======================================
        flux_1_m = -mu_m*grad1_u3_m
        flux_2_m = -mu_m*grad2_u3_m
        flux_3_m = -mu_m*grad3_u3_m

        flux_1_p = -mu_p*grad1_u3_p
        flux_2_p = -mu_p*grad2_u3_p
        flux_3_p = -mu_p*grad3_u3_p


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_average('grid_displacement3','Diffusion',    &
                                                flux_1_m, flux_2_m, flux_3_m,       &
                                                flux_1_p, flux_2_p, flux_3_p)


    end subroutine compute
    !**************************************************************************************************




end module MMD_boundary_operator
