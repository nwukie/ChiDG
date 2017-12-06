module graddemo_P_boundary_operator
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
    type, extends(operator_t), public :: graddemo_P_boundary_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type graddemo_P_boundary_operator_t
    !********************************************************************************



contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(graddemo_P_boundary_operator_t),   intent(inout) :: self
        
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
        class(graddemo_P_boundary_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop



        type(AD_D), allocatable, dimension(:)   ::  &
            grad1_pbc_m, grad2_pbc_m, grad3_pbc_m,  &
            grad1_pbc_p, grad2_pbc_p, grad3_pbc_p,  &
            grad1_p_m,   grad2_p_m,   grad3_p_m,    &
            grad1_p_p,   grad2_p_p,   grad3_p_p,    &
            flux_1_m, flux_2_m, flux_3_m,           &
            flux_1_p, flux_2_p, flux_3_p

    

        !
        ! Interpolate solution to quadrature nodes
        !
        grad1_pbc_m = worker%get_field('Pressure', 'grad1', 'face interior')
        grad2_pbc_m = worker%get_field('Pressure', 'grad2', 'face interior')
        grad3_pbc_m = worker%get_field('Pressure', 'grad3', 'face interior')

        grad1_pbc_p = worker%get_field('Pressure', 'grad1', 'face exterior')
        grad2_pbc_p = worker%get_field('Pressure', 'grad2', 'face exterior')
        grad3_pbc_p = worker%get_field('Pressure', 'grad3', 'face exterior')


!        grad1_p_m =
!        grad2_p_m =
!        grad3_p_m =

!        grad1_p_p =
!        grad2_p_p =
!        grad3_p_p =





        flux_1_m = grad1_pbc_m - grad1_p_m
        flux_2_m = grad2_pbc_m - grad2_p_m
        flux_3_m = grad3_pbc_m - grad3_p_m

        flux_1_p = grad1_pbc_p - grad1_p_p
        flux_2_p = grad2_pbc_p - grad2_p_p
        flux_3_p = grad3_pbc_p - grad3_p_p


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_average('Pressure','Diffusion',          &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)


    end subroutine compute
    !**************************************************************************************************




end module graddemo_P_boundary_operator
