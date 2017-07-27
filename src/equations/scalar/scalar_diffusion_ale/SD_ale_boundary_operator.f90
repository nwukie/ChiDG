module SD_ale_boundary_operator
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
    type, extends(operator_t), public :: SD_ale_boundary_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type SD_ale_boundary_operator_t
    !********************************************************************************



contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(SD_ale_boundary_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Scalar Diffusion ALE Boundary Average Operator")

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
        class(SD_ale_boundary_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop



        type(AD_D), allocatable, dimension(:)   ::  &
            flux_1, flux_2, flux_3,                 &
            flux_m, flux_p, integrand,        &
            mu_m, mu_p

        real(rk),   allocatable, dimension(:)   :: &
            norm_1, norm_2, norm_3

        type(AD_D), allocatable, dimension(:,:)   ::  &
            gradu_m, gradu_p, &
            flux_ref



        !
        ! Interpolate solution to quadrature nodes
        !
        gradu_m = worker%get_primary_field_grad_ale_face('u', 'gradient + lift','face interior')
        gradu_p = worker%get_primary_field_grad_ale_face('u', 'gradient + lift','face exterior')


        norm_1 = worker%normal(1)
        norm_2 = worker%normal(2)
        norm_3 = worker%normal(3)


        !
        ! Compute scalar coefficient
        !
        mu_m = worker%get_model_field_face('Scalar Diffusion Coefficient', 'value', 'face interior')
        mu_p = worker%get_model_field_face('Scalar Diffusion Coefficient', 'value', 'face exterior')


        flux_m = -mu_m*gradu_m(:,1)
        flux_p = -mu_p*gradu_p(:,1)
        flux_1 = HALF*(flux_m + flux_p)

        flux_m = -mu_m*gradu_m(:,2)
        flux_p = -mu_p*gradu_p(:,2)
        flux_2 = HALF*(flux_m + flux_p)

        flux_m = -mu_m*gradu_m(:,3)
        flux_p = -mu_p*gradu_p(:,3)
        flux_3 = HALF*(flux_m + flux_p)


        flux_ref = worker%post_process_boundary_diffusive_flux_ale(flux_1, flux_2, flux_3, 'face interior')
        !
        ! Compute boundary average flux
        !
        integrand = flux_ref(:,1)*norm_1 + flux_ref(:,2)*norm_2 + flux_ref(:,3)*norm_3


        !
        ! Integrate flux
        !
        call worker%integrate_boundary('u',integrand)


    end subroutine compute
    !**************************************************************************************************




end module SD_ale_boundary_operator
