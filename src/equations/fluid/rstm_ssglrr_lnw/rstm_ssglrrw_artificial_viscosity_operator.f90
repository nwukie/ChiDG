module rstm_ssglrrw_artificial_viscosity_operator
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF

    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    use mod_rstm_ssglrrw, only: rstm_ssglrrw_avc
    use ieee_arithmetic
    implicit none

    private

    
    !> Volume flux for Fluid laplacian_av Terms.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/23/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: rstm_ssglrrw_artificial_viscosity_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_artificial_viscosity_operator_t
    !******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/23/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rstm_ssglrrw_artificial_viscosity_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name('RSTMSSGLRRW Artificial Viscosity Volume Operator')

        ! Set operator type
        call self%set_operator_type('Volume Diffusive Operator')

        ! Set operator equations
        call self%add_primary_field('Density * Reynolds-11'   )
        call self%add_primary_field('Density * Reynolds-22'   )
        call self%add_primary_field('Density * Reynolds-33'   )
        call self%add_primary_field('Density * Reynolds-12'   )
        call self%add_primary_field('Density * Reynolds-13'   )
        call self%add_primary_field('Density * Reynolds-23'   )

    end subroutine init
    !********************************************************************************



    !> Volume flux routine for Fluid laplacian_av Terms.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/23/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(rstm_ssglrrw_artificial_viscosity_operator_t), intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:) ::    &
            mu_neg, grad1, grad2, grad3, &
            flux_1, flux_2, flux_3

        real(rk),   allocatable, dimension(:)   :: r

        integer(ik) :: p
        real(rk)    :: h(3), sigma
        real(rk),   allocatable :: h_smooth(:, :)


       
        !
        ! Get Model fields:
        !   Second Coefficient of Viscosity
        !   Thermal Conductivity
        !

        h_smooth = worker%h_smooth('element')

        
        h = worker%element_size('interior')
        p = worker%solution_order('interior')
        h = h/real(p+1,rk)
        h_smooth = h_smooth/real(p+1,rk)
        sigma = 0.5_rk

        mu_neg = worker%get_field('RSTM AV-11'   , 'value', 'element')
        grad1 = worker%get_field('Density * Reynolds-11'   , 'grad1', 'element')
        grad2 = worker%get_field('Density * Reynolds-11'   , 'grad2', 'element')
        grad3 = worker%get_field('Density * Reynolds-11'   , 'grad3', 'element')
        !flux_1 = -(h(1)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad1
        !flux_2 = -(h(2)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad2
        !flux_3 = -(h(3)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad3
        
        flux_1 = -(h_smooth(:,1)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad1
        flux_2 = -(h_smooth(:,2)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad2
        flux_3 = -(h_smooth(:,3)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad3
        
        call worker%integrate_volume_flux('Density * Reynolds-11','Diffusion',flux_1,flux_2,flux_3)

        mu_neg = worker%get_field('RSTM AV-22'   , 'value', 'element')
        grad1 = worker%get_field('Density * Reynolds-22'   , 'grad1', 'element')
        grad2 = worker%get_field('Density * Reynolds-22'   , 'grad2', 'element')
        grad3 = worker%get_field('Density * Reynolds-22'   , 'grad3', 'element')
        !flux_1 = -(h(1)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad1
        !flux_2 = -(h(2)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad2
        !flux_3 = -(h(3)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad3
        
        flux_1 = -(h_smooth(:,1)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad1
        flux_2 = -(h_smooth(:,2)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad2
        flux_3 = -(h_smooth(:,3)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad3
        
        call worker%integrate_volume_flux('Density * Reynolds-22','Diffusion',flux_1,flux_2,flux_3)

        mu_neg = worker%get_field('RSTM AV-33'   , 'value', 'element')
        grad1 = worker%get_field('Density * Reynolds-33'   , 'grad1', 'element')
        grad2 = worker%get_field('Density * Reynolds-33'   , 'grad2', 'element')
        grad3 = worker%get_field('Density * Reynolds-33'   , 'grad3', 'element')
        !flux_1 = -(h(1)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad1
        !flux_2 = -(h(2)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad2
        !flux_3 = -(h(3)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad3
        
        flux_1 = -(h_smooth(:,1)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad1
        flux_2 = -(h_smooth(:,2)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad2
        flux_3 = -(h_smooth(:,3)**TWO*rstm_ssglrrw_avc+mu_neg*sigma)*grad3
        
        call worker%integrate_volume_flux('Density * Reynolds-33','Diffusion',flux_1,flux_2,flux_3)

        grad1 = worker%get_field('Density * Reynolds-12'   , 'grad1', 'element')
        grad2 = worker%get_field('Density * Reynolds-12'   , 'grad2', 'element')
        grad3 = worker%get_field('Density * Reynolds-12'   , 'grad3', 'element')
        !flux_1 = -(h(1)**TWO*rstm_ssglrrw_avc)*grad1
        !flux_2 = -(h(2)**TWO*rstm_ssglrrw_avc)*grad2
        !flux_3 = -(h(3)**TWO*rstm_ssglrrw_avc)*grad3
        
        
        flux_1 = -(h_smooth(:,1)**TWO*rstm_ssglrrw_avc)*grad1
        flux_2 = -(h_smooth(:,2)**TWO*rstm_ssglrrw_avc)*grad2
        flux_3 = -(h_smooth(:,3)**TWO*rstm_ssglrrw_avc)*grad3
        call worker%integrate_volume_flux('Density * Reynolds-12','Diffusion',flux_1,flux_2,flux_3)

        grad1 = worker%get_field('Density * Reynolds-13'   , 'grad1', 'element')
        grad2 = worker%get_field('Density * Reynolds-13'   , 'grad2', 'element')
        grad3 = worker%get_field('Density * Reynolds-13'   , 'grad3', 'element')
        !flux_1 = -(h(1)**TWO*rstm_ssglrrw_avc)*grad1
        !flux_2 = -(h(2)**TWO*rstm_ssglrrw_avc)*grad2
        !flux_3 = -(h(3)**TWO*rstm_ssglrrw_avc)*grad3
        
        
        flux_1 = -(h_smooth(:,1)**TWO*rstm_ssglrrw_avc)*grad1
        flux_2 = -(h_smooth(:,2)**TWO*rstm_ssglrrw_avc)*grad2
        flux_3 = -(h_smooth(:,3)**TWO*rstm_ssglrrw_avc)*grad3
        call worker%integrate_volume_flux('Density * Reynolds-13','Diffusion',flux_1,flux_2,flux_3)

        grad1 = worker%get_field('Density * Reynolds-23'   , 'grad1', 'element')
        grad2 = worker%get_field('Density * Reynolds-23'   , 'grad2', 'element')
        grad3 = worker%get_field('Density * Reynolds-23'   , 'grad3', 'element')
        !flux_1 = -(h(1)**TWO*rstm_ssglrrw_avc)*grad1
        !flux_2 = -(h(2)**TWO*rstm_ssglrrw_avc)*grad2
        !flux_3 = -(h(3)**TWO*rstm_ssglrrw_avc)*grad3
        
        
        flux_1 = -(h_smooth(:,1)**TWO*rstm_ssglrrw_avc)*grad1
        flux_2 = -(h_smooth(:,2)**TWO*rstm_ssglrrw_avc)*grad2
        flux_3 = -(h_smooth(:,3)**TWO*rstm_ssglrrw_avc)*grad3
        call worker%integrate_volume_flux('Density * Reynolds-23','Diffusion',flux_1,flux_2,flux_3)





    end subroutine compute
    !*********************************************************************************************************






end module rstm_ssglrrw_artificial_viscosity_operator
