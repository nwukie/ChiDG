module rstm_ssglrrw_artificial_viscosity_boundary_average_operator
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, HALF

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    use mod_rstm_ssglrrw, only: rstm_ssglrrw_avc
    use ieee_arithmetic
    implicit none

    private



    !> Implementation of the Euler boundary average flux
    !!
    !!  - At a boundary interface, the solution states Q- and Q+ exists on opposite 
    !!    sides of the boundary. The average flux is computed as Favg = 1/2(F(Q-) + F(Q+))
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: rstm_ssglrrw_artificial_viscosity_boundary_average_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rstm_ssglrrw_artificial_viscosity_boundary_average_operator_t
    !********************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rstm_ssglrrw_artificial_viscosity_boundary_average_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('RSTMSSGLRRW Artificial Viscosity Boundary Average Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('Boundary Diffusive Flux')

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



    !>  Boundary Flux routine for Euler
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!-------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(rstm_ssglrrw_artificial_viscosity_boundary_average_operator_t),   intent(inout)   :: self
        type(chidg_worker_t),                               intent(inout)   :: worker
        class(properties_t),                                intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::                        &
            grad1_m, grad2_m, grad3_m, pgrad1_m, pgrad2_m, pgrad3_m,    &
            grad1_p, grad2_p, grad3_p, pgrad1_p, pgrad2_p, pgrad3_p,    &
            mu_neg_m, mu_neg_p, sigma_k_m, sigma_k_p,                                                 &
            flux_1_m, flux_2_m, flux_3_m,                               &
            flux_1_p, flux_2_p, flux_3_p
        
        integer(ik) :: p
        real(rk) :: sigma, h_m(3), h_p(3)
        real(rk),   allocatable :: h_s_m(:,:), h_s_p(:,:)
    
        sigma = 0.5_rk


        h_s_m = worker%h_smooth('face interior')
        h_s_p = worker%h_smooth('face exterior')
        h_m = worker%element_size('interior')
        h_p = worker%element_size('exterior')
        p = worker%solution_order('interior')
        h_m = h_m/real(p+1,rk)
        h_p = h_p/real(p+1,rk)
        h_s_m = h_s_m/real(p+1,rk)
        h_s_p = h_s_p/real(p+1,rk)

        mu_neg_m = worker%get_field('RSTM AV-11', 'value', 'face interior')
        mu_neg_p = worker%get_field('RSTM AV-11', 'value', 'face exterior')
        !----------------------------------
        ! 11            
        !----------------------------------
        grad1_m = worker%get_field('Density * Reynolds-11'   , 'grad1', 'face interior')
        grad2_m = worker%get_field('Density * Reynolds-11'   , 'grad2', 'face interior')
        grad3_m = worker%get_field('Density * Reynolds-11'   , 'grad3', 'face interior')
        flux_1_m = -(h_s_m(:,1)**TWO*rstm_ssglrrw_avc+sigma*mu_neg_m)*grad1_m
        flux_2_m = -(h_s_m(:,2)**TWO*rstm_ssglrrw_avc+sigma*mu_neg_m)*grad2_m
        flux_3_m = -(h_s_m(:,3)**TWO*rstm_ssglrrw_avc+sigma*mu_neg_m)*grad3_m

        grad1_p = worker%get_field('Density * Reynolds-11'   , 'grad1', 'face exterior')
        grad2_p = worker%get_field('Density * Reynolds-11'   , 'grad2', 'face exterior')
        grad3_p = worker%get_field('Density * Reynolds-11'   , 'grad3', 'face exterior')
        flux_1_p = -(h_s_p(:,1)**TWO*rstm_ssglrrw_avc + sigma*mu_neg_p)*grad1_p
        flux_2_p = -(h_s_p(:,2)**TWO*rstm_ssglrrw_avc + sigma*mu_neg_p)*grad2_p
        flux_3_p = -(h_s_p(:,3)**TWO*rstm_ssglrrw_avc + sigma*mu_neg_p)*grad3_p


        call worker%integrate_boundary_average('Density * Reynolds-11','Diffusion',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)

        mu_neg_m = worker%get_field('RSTM AV-22', 'value', 'face interior')
        mu_neg_p = worker%get_field('RSTM AV-22', 'value', 'face exterior')
        !----------------------------------
        !  22          
        !----------------------------------
        grad1_m = worker%get_field('Density * Reynolds-22'   , 'grad1', 'face interior')
        grad2_m = worker%get_field('Density * Reynolds-22'   , 'grad2', 'face interior')
        grad3_m = worker%get_field('Density * Reynolds-22'   , 'grad3', 'face interior')
        flux_1_m = -(h_s_m(:,1)**TWO*rstm_ssglrrw_avc+sigma*mu_neg_m)*grad1_m
        flux_2_m = -(h_s_m(:,2)**TWO*rstm_ssglrrw_avc+sigma*mu_neg_m)*grad2_m
        flux_3_m = -(h_s_m(:,3)**TWO*rstm_ssglrrw_avc+sigma*mu_neg_m)*grad3_m

        grad1_p = worker%get_field('Density * Reynolds-22'   , 'grad1', 'face exterior')
        grad2_p = worker%get_field('Density * Reynolds-22'   , 'grad2', 'face exterior')
        grad3_p = worker%get_field('Density * Reynolds-22'   , 'grad3', 'face exterior')
        flux_1_p = -(h_s_p(:,1)**TWO*rstm_ssglrrw_avc + sigma*mu_neg_p)*grad1_p
        flux_2_p = -(h_s_p(:,2)**TWO*rstm_ssglrrw_avc + sigma*mu_neg_p)*grad2_p
        flux_3_p = -(h_s_p(:,3)**TWO*rstm_ssglrrw_avc + sigma*mu_neg_p)*grad3_p


        call worker%integrate_boundary_average('Density * Reynolds-22','Diffusion',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)


        mu_neg_m = worker%get_field('RSTM AV-33', 'value', 'face interior')
        mu_neg_p = worker%get_field('RSTM AV-33', 'value', 'face exterior')
        !----------------------------------
        !           33 
        !----------------------------------
        grad1_m = worker%get_field('Density * Reynolds-33'   , 'grad1', 'face interior')
        grad2_m = worker%get_field('Density * Reynolds-33'   , 'grad2', 'face interior')
        grad3_m = worker%get_field('Density * Reynolds-33'   , 'grad3', 'face interior')
        flux_1_m = -(h_s_m(:,1)**TWO*rstm_ssglrrw_avc+sigma*mu_neg_m)*grad1_m
        flux_2_m = -(h_s_m(:,2)**TWO*rstm_ssglrrw_avc+sigma*mu_neg_m)*grad2_m
        flux_3_m = -(h_s_m(:,3)**TWO*rstm_ssglrrw_avc+sigma*mu_neg_m)*grad3_m

        grad1_p = worker%get_field('Density * Reynolds-33'   , 'grad1', 'face exterior')
        grad2_p = worker%get_field('Density * Reynolds-33'   , 'grad2', 'face exterior')
        grad3_p = worker%get_field('Density * Reynolds-33'   , 'grad3', 'face exterior')
        flux_1_p = -(h_s_p(:,1)**TWO*rstm_ssglrrw_avc + sigma*mu_neg_p)*grad1_p
        flux_2_p = -(h_s_p(:,2)**TWO*rstm_ssglrrw_avc + sigma*mu_neg_p)*grad2_p
        flux_3_p = -(h_s_p(:,3)**TWO*rstm_ssglrrw_avc + sigma*mu_neg_p)*grad3_p


        call worker%integrate_boundary_average('Density * Reynolds-33','Diffusion',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)

        !----------------------------------
        !           12 
        !----------------------------------
        grad1_m = worker%get_field('Density * Reynolds-12'   , 'grad1', 'face interior')
        grad2_m = worker%get_field('Density * Reynolds-12'   , 'grad2', 'face interior')
        grad3_m = worker%get_field('Density * Reynolds-12'   , 'grad3', 'face interior')
        flux_1_m = -(h_s_m(:,1)**TWO*rstm_ssglrrw_avc)*grad1_m
        flux_2_m = -(h_s_m(:,2)**TWO*rstm_ssglrrw_avc)*grad2_m
        flux_3_m = -(h_s_m(:,3)**TWO*rstm_ssglrrw_avc)*grad3_m

        grad1_p = worker%get_field('Density * Reynolds-12'   , 'grad1', 'face exterior')
        grad2_p = worker%get_field('Density * Reynolds-12'   , 'grad2', 'face exterior')
        grad3_p = worker%get_field('Density * Reynolds-12'   , 'grad3', 'face exterior')
        flux_1_p = -(h_s_p(:,1)**TWO*rstm_ssglrrw_avc)*grad1_p
        flux_2_p = -(h_s_p(:,2)**TWO*rstm_ssglrrw_avc)*grad2_p
        flux_3_p = -(h_s_p(:,3)**TWO*rstm_ssglrrw_avc)*grad3_p


        call worker%integrate_boundary_average('Density * Reynolds-12','Diffusion',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)

        !----------------------------------
        !           13 
        !----------------------------------
        grad1_m = worker%get_field('Density * Reynolds-13'   , 'grad1', 'face interior')
        grad2_m = worker%get_field('Density * Reynolds-13'   , 'grad2', 'face interior')
        grad3_m = worker%get_field('Density * Reynolds-13'   , 'grad3', 'face interior')
        flux_1_m = -(h_s_m(:,1)**TWO*rstm_ssglrrw_avc)*grad1_m
        flux_2_m = -(h_s_m(:,2)**TWO*rstm_ssglrrw_avc)*grad2_m
        flux_3_m = -(h_s_m(:,3)**TWO*rstm_ssglrrw_avc)*grad3_m

        grad1_p = worker%get_field('Density * Reynolds-13'   , 'grad1', 'face exterior')
        grad2_p = worker%get_field('Density * Reynolds-13'   , 'grad2', 'face exterior')
        grad3_p = worker%get_field('Density * Reynolds-13'   , 'grad3', 'face exterior')
        flux_1_p = -(h_s_p(:,1)**TWO*rstm_ssglrrw_avc)*grad1_p
        flux_2_p = -(h_s_p(:,2)**TWO*rstm_ssglrrw_avc)*grad2_p
        flux_3_p = -(h_s_p(:,3)**TWO*rstm_ssglrrw_avc)*grad3_p


        call worker%integrate_boundary_average('Density * Reynolds-13','Diffusion',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)

        !----------------------------------
        !           23 
        !----------------------------------
        grad1_m = worker%get_field('Density * Reynolds-23'   , 'grad1', 'face interior')
        grad2_m = worker%get_field('Density * Reynolds-23'   , 'grad2', 'face interior')
        grad3_m = worker%get_field('Density * Reynolds-23'   , 'grad3', 'face interior')
        flux_1_m = -(h_s_m(:,1)**TWO*rstm_ssglrrw_avc)*grad1_m
        flux_2_m = -(h_s_m(:,2)**TWO*rstm_ssglrrw_avc)*grad2_m
        flux_3_m = -(h_s_m(:,3)**TWO*rstm_ssglrrw_avc)*grad3_m

        grad1_p = worker%get_field('Density * Reynolds-23'   , 'grad1', 'face exterior')
        grad2_p = worker%get_field('Density * Reynolds-23'   , 'grad2', 'face exterior')
        grad3_p = worker%get_field('Density * Reynolds-23'   , 'grad3', 'face exterior')
        flux_1_p = -(h_s_p(:,1)**TWO*rstm_ssglrrw_avc)*grad1_p
        flux_2_p = -(h_s_p(:,2)**TWO*rstm_ssglrrw_avc)*grad2_p
        flux_3_p = -(h_s_p(:,3)**TWO*rstm_ssglrrw_avc)*grad3_p


        call worker%integrate_boundary_average('Density * Reynolds-23','Diffusion',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)




    end subroutine compute
    !*********************************************************************************************************












end module rstm_ssglrrw_artificial_viscosity_boundary_average_operator
