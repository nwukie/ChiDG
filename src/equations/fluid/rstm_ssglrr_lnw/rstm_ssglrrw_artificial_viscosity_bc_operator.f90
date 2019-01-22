module rstm_ssglrrw_artificial_viscosity_bc_operator
#include <messenger.h>
    use mod_kinds,          only: ik, rk
    use mod_constants,      only: HALF, ONE, TWO
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    use mod_rstm_ssglrrw, only: rstm_ssglrrw_avc
    use ieee_arithmetic
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !------------------------------------------------------------------------------------
    type, public, extends(operator_t)   :: rstm_ssglrrw_artificial_viscosity_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_artificial_viscosity_bc_operator_t
    !*************************************************************************************






contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rstm_ssglrrw_artificial_viscosity_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('RSTMSSGLRRW Artificial Viscosity BC Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('BC Diffusive Flux')

        !
        ! Set operator equations
        !
        call self%add_primary_field('Density * Reynolds-11'   )
        call self%add_primary_field('Density * Reynolds-22'   )
        call self%add_primary_field('Density * Reynolds-33'   )
        call self%add_primary_field('Density * Reynolds-12'   )
        call self%add_primary_field('Density * Reynolds-13'   )
        call self%add_primary_field('Density * Reynolds-23'   )
    end subroutine init
    !********************************************************************************





    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      mesh    Mesh data containing boundarys and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!  @param[in]      face    face_info_t containing indices on location and misc information
    !!  @param[in]      flux    function_into_t containing info on the function being computed
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(rstm_ssglrrw_artificial_viscosity_bc_operator_t), intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::    &
            mu_neg, grad1, grad2, grad3, &
            flux_1, flux_2, flux_3


        integer(ik) :: p
        real(rk) :: sigma, h(3)

               

        h = worker%element_size('interior')
        p = worker%solution_order('interior')
        h = h/real(p+1,rk)
        sigma = 0.5_rk

        mu_neg = worker%get_field('RSTM AV-11'   , 'value', 'boundary')
        grad1 = worker%get_field('Density * Reynolds-11'   , 'grad1', 'boundary')
        grad2 = worker%get_field('Density * Reynolds-11'   , 'grad2', 'boundary')
        grad3 = worker%get_field('Density * Reynolds-11'   , 'grad3', 'boundary')
        flux_1 = -(h(1)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad1
        flux_2 = -(h(2)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad2
        flux_3 = -(h(3)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad3
        
        
        call worker%integrate_boundary_condition('Density * Reynolds-11','Diffusion',flux_1,flux_2,flux_3)

        mu_neg = worker%get_field('RSTM AV-22'   , 'value', 'boundary')
        grad1 = worker%get_field('Density * Reynolds-22'   , 'grad1', 'boundary')
        grad2 = worker%get_field('Density * Reynolds-22'   , 'grad2', 'boundary')
        grad3 = worker%get_field('Density * Reynolds-22'   , 'grad3', 'boundary')
        flux_1 = -(h(1)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad1
        flux_2 = -(h(2)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad2
        flux_3 = -(h(3)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad3
        
        
        call worker%integrate_boundary_condition('Density * Reynolds-22','Diffusion',flux_1,flux_2,flux_3)

        mu_neg = worker%get_field('RSTM AV-33'   , 'value', 'boundary')
        grad1 = worker%get_field('Density * Reynolds-33'   , 'grad1', 'boundary')
        grad2 = worker%get_field('Density * Reynolds-33'   , 'grad2', 'boundary')
        grad3 = worker%get_field('Density * Reynolds-33'   , 'grad3', 'boundary')
        flux_1 = -(h(1)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad1
        flux_2 = -(h(2)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad2
        flux_3 = -(h(3)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad3
        
        
        call worker%integrate_boundary_condition('Density * Reynolds-33','Diffusion',flux_1,flux_2,flux_3)

        grad1 = worker%get_field('Density * Reynolds-12'   , 'grad1', 'boundary')
        grad2 = worker%get_field('Density * Reynolds-12'   , 'grad2', 'boundary')
        grad3 = worker%get_field('Density * Reynolds-12'   , 'grad3', 'boundary')
        flux_1 = -(h(1)**TWO*rstm_ssglrrw_avc)*grad1
        flux_2 = -(h(2)**TWO*rstm_ssglrrw_avc)*grad2
        flux_3 = -(h(3)**TWO*rstm_ssglrrw_avc)*grad3
        
        
        call worker%integrate_boundary_condition('Density * Reynolds-12','Diffusion',flux_1,flux_2,flux_3)

        grad1 = worker%get_field('Density * Reynolds-13'   , 'grad1', 'boundary')
        grad2 = worker%get_field('Density * Reynolds-13'   , 'grad2', 'boundary')
        grad3 = worker%get_field('Density * Reynolds-13'   , 'grad3', 'boundary')
        flux_1 = -(h(1)**TWO*rstm_ssglrrw_avc)*grad1
        flux_2 = -(h(2)**TWO*rstm_ssglrrw_avc)*grad2
        flux_3 = -(h(3)**TWO*rstm_ssglrrw_avc)*grad3
        
        
        call worker%integrate_boundary_condition('Density * Reynolds-13','Diffusion',flux_1,flux_2,flux_3)

        grad1 = worker%get_field('Density * Reynolds-23'   , 'grad1', 'boundary')
        grad2 = worker%get_field('Density * Reynolds-23'   , 'grad2', 'boundary')
        grad3 = worker%get_field('Density * Reynolds-23'   , 'grad3', 'boundary')
        flux_1 = -(h(1)**TWO*rstm_ssglrrw_avc)*grad1
        flux_2 = -(h(2)**TWO*rstm_ssglrrw_avc)*grad2
        flux_3 = -(h(3)**TWO*rstm_ssglrrw_avc)*grad3
        
        
        call worker%integrate_boundary_condition('Density * Reynolds-23','Diffusion',flux_1,flux_2,flux_3)





    end subroutine compute
    !**********************************************************************************************









end module rstm_ssglrrw_artificial_viscosity_bc_operator
