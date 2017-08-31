module SD_boundary_operator
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
    type, extends(operator_t), public :: SD_boundary_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type SD_boundary_operator_t
    !********************************************************************************



contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(SD_boundary_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Scalar Diffusion Boundary Average Operator")

        !
        ! Set operator type
        !
        call self%set_operator_type("Boundary Diffusive Operator")

        !
        ! Set operator equations
        !
        call self%add_primary_field("u")

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
        class(SD_boundary_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop



        type(AD_D), allocatable, dimension(:)   ::  &
            grad1_u_m, grad2_u_m, grad3_u_m,        &
            grad1_u_p, grad2_u_p, grad3_u_p,        &
            flux_1_m, flux_2_m, flux_3_m,           &
            flux_1_p, flux_2_p, flux_3_p,           &
            mu_m, mu_p


        !
        ! Interpolate solution to quadrature nodes
        !
        grad1_u_m = worker%get_field('u', 'grad1', 'face interior')
        grad2_u_m = worker%get_field('u', 'grad2', 'face interior')
        grad3_u_m = worker%get_field('u', 'grad3', 'face interior')

        grad1_u_p = worker%get_field('u', 'grad1', 'face exterior')
        grad2_u_p = worker%get_field('u', 'grad2', 'face exterior')
        grad3_u_p = worker%get_field('u', 'grad3', 'face exterior')


        !
        ! Compute scalar coefficient
        !
        mu_m = worker%get_field('Scalar Diffusion Coefficient', 'value', 'face interior')
        mu_p = worker%get_field('Scalar Diffusion Coefficient', 'value', 'face exterior')


        flux_1_m = -mu_m*grad1_u_m
        flux_2_m = -mu_m*grad2_u_m
        flux_3_m = -mu_m*grad3_u_m

        flux_1_p = -mu_p*grad1_u_p
        flux_2_p = -mu_p*grad2_u_p
        flux_3_p = -mu_p*grad3_u_p


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_average('u','Diffusion',                 &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)


    end subroutine compute
    !**************************************************************************************************




end module SD_boundary_operator
