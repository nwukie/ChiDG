module sst_artificial_viscosity_boundary_average_operator
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, HALF

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    use mod_sst, only: sst_avc
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
    type, extends(operator_t), public :: sst_artificial_viscosity_boundary_average_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type sst_artificial_viscosity_boundary_average_operator_t
    !********************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(sst_artificial_viscosity_boundary_average_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('SST Artificial Viscosity Boundary Average Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('Boundary Diffusive Flux')

        !
        ! Set operator equations
        !
        call self%add_primary_field('Density * k'   )
    end subroutine init
    !********************************************************************************



    !>  Boundary Flux routine for Euler
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!-------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(sst_artificial_viscosity_boundary_average_operator_t),   intent(inout)   :: self
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
        real(rk) :: avc, h_m(3), h_p(3)
    
        sigma_k_m = worker%get_field('SST sigma_k', 'value', 'face interior')
        sigma_k_p = worker%get_field('SST sigma_k', 'value', 'face exterior')

        mu_neg_m = worker%get_field('SST Artificial Viscosity k-neg', 'value', 'face interior')
        mu_neg_p = worker%get_field('SST Artificial Viscosity k-neg', 'value', 'face exterior')

        h_m = worker%element_size('interior')
        h_p = worker%element_size('exterior')
        p = worker%solution_order('interior')
        h_m = h_m/real(p+1,rk)
        h_p = h_p/real(p+1,rk)

        !----------------------------------
        !            mass flux
        !----------------------------------
        grad1_m = worker%get_field('Density * k'   , 'grad1', 'face interior')
        grad2_m = worker%get_field('Density * k'   , 'grad2', 'face interior')
        grad3_m = worker%get_field('Density * k'   , 'grad3', 'face interior')
        flux_1_m = -(h_m(1)**TWO*sst_avc+sigma_k_m*mu_neg_m)*grad1_m
        flux_2_m = -(h_m(2)**TWO*sst_avc+sigma_k_m*mu_neg_m)*grad2_m
        flux_3_m = -(h_m(3)**TWO*sst_avc+sigma_k_m*mu_neg_m)*grad3_m

        grad1_p = worker%get_field('Density * k'   , 'grad1', 'face exterior')
        grad2_p = worker%get_field('Density * k'   , 'grad2', 'face exterior')
        grad3_p = worker%get_field('Density * k'   , 'grad3', 'face exterior')
        flux_1_p = -(h_p(1)**TWO*sst_avc + sigma_k_p*mu_neg_p)*grad1_p
        flux_2_p = -(h_p(2)**TWO*sst_avc + sigma_k_p*mu_neg_p)*grad2_p
        flux_3_p = -(h_p(3)**TWO*sst_avc + sigma_k_p*mu_neg_p)*grad3_p


        call worker%integrate_boundary_average('Density * k','Diffusion',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)


    end subroutine compute
    !*********************************************************************************************************












end module sst_artificial_viscosity_boundary_average_operator
