module spalart_allmaras_laxfriedrichs
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF
    use mod_fluid,              only: omega
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    type, extends(operator_t), public :: spalart_allmaras_laxfriedrichs_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type spalart_allmaras_laxfriedrichs_operator_t
    !*****************************************************************************************

contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init(self)
        class(spalart_allmaras_laxfriedrichs_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('Spalart-Allmaras LaxFriedrichs Operator')

        ! Set operator type
        call self%set_operator_type('Boundary Advective Operator')

        ! Set operator equations
        call self%add_primary_field('Density * NuTilde')

    end subroutine init
    !*****************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/9/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(spalart_allmaras_laxfriedrichs_operator_t), intent(inout)   :: self
        type(chidg_worker_t),                             intent(inout)   :: worker
        class(properties_t),                              intent(inout)   :: prop


        type(AD_D), dimension(:), allocatable   ::          &
            density_nutilde_m, density_nutilde_p,           &
            u_a_m, v_a_m, w_a_m, T_m, un_m, c_m, diss_m,    &
            u_a_p, v_a_p, w_a_p, T_p, un_p, c_p, diss_p,    &
            flux_avg_1, flux_avg_2, flux_avg_3, diff,       &
            flux_1, flux_2, flux_3, integrand

        real(rk),   dimension(:), allocatable   ::  &
            norm_1,  norm_2,  norm_3,               &
            unorm_1, unorm_2, unorm_3, area



        !
        ! Interpolate solution to quadrature nodes
        !
        density_nutilde_m = worker%get_primary_field_face('Density * NuTilde', 'value', 'face interior')
        density_nutilde_p = worker%get_primary_field_face('Density * NuTilde', 'value', 'face exterior')


        
        !
        ! Get fluid advection velocity
        !
        u_a_m = worker%get_model_field_face('Advection Velocity-1', 'value', 'face interior')
        v_a_m = worker%get_model_field_face('Advection Velocity-2', 'value', 'face interior')
        w_a_m = worker%get_model_field_face('Advection Velocity-3', 'value', 'face interior')

        u_a_p = worker%get_model_field_face('Advection Velocity-1', 'value', 'face exterior')
        v_a_p = worker%get_model_field_face('Advection Velocity-2', 'value', 'face exterior')
        w_a_p = worker%get_model_field_face('Advection Velocity-3', 'value', 'face exterior')


        !
        ! Get normal vector
        !
        norm_1  = worker%normal(1)
        norm_2  = worker%normal(2)
        norm_3  = worker%normal(3)
        unorm_1 = worker%unit_normal(1)
        unorm_2 = worker%unit_normal(2)
        unorm_3 = worker%unit_normal(3)


        !
        ! Compute average flux
        ! 
        flux_avg_1 = HALF*(density_nutilde_m*u_a_m  +  density_nutilde_p*u_a_p)
        flux_avg_2 = HALF*(density_nutilde_m*v_a_m  +  density_nutilde_p*v_a_p)
        flux_avg_3 = HALF*(density_nutilde_m*w_a_m  +  density_nutilde_p*w_a_p)


        !
        ! Compute maximum wave speed
        !
        un_m = u_a_m*unorm_1 + v_a_m*unorm_2 + w_a_m*unorm_3
        un_p = u_a_p*unorm_1 + v_a_p*unorm_2 + w_a_p*unorm_3


        T_m = worker%get_model_field_face('Temperature','value','face interior')
        T_p = worker%get_model_field_face('Temperature','value','face exterior')
        c_m = sqrt(1.4_rk * 287.15_rk * T_m)
        c_p = sqrt(1.4_rk * 287.15_rk * T_p)

        diss_m = abs(un_m) + c_m
        diss_p = abs(un_p) + c_p


        !
        ! Compute Lax-Friedrichs upwind flux
        !
        diff   = (density_nutilde_m - density_nutilde_p)
        flux_1 = flux_avg_1 + max(abs(diss_m),abs(diss_p))*HALF*diff
        flux_2 = flux_avg_2 + max(abs(diss_m),abs(diss_p))*HALF*diff
        flux_3 = flux_avg_3 + max(abs(diss_m),abs(diss_p))*HALF*diff


        integrand = flux_avg_1*norm_1 + flux_avg_2*norm_2 + flux_avg_3*norm_3

        integrand = integrand + max(abs(diss_m),abs(diss_p))*HALF*diff*norm_1*unorm_1
        integrand = integrand + max(abs(diss_m),abs(diss_p))*HALF*diff*norm_2*unorm_2
        integrand = integrand + max(abs(diss_m),abs(diss_p))*HALF*diff*norm_3*unorm_3


!        integrand = flux_avg_1*norm_1 + flux_avg_2*norm_2 + flux_avg_3*norm_3
!
!        diff   = (density_nutilde_p - density_nutilde_m)
!        area = sqrt(norm_1**TWO + norm_2**TWO + norm_3**TWO)
!        integrand = integrand + max(abs(diss_m),abs(diss_p))*HALF*diff*area
!        integrand = integrand + max(abs(diss_m),abs(diss_p))*HALF*diff*area
!        integrand = integrand + max(abs(diss_m),abs(diss_p))*HALF*diff*area

        !
        ! Integrate flux
        !
        call worker%integrate_boundary('Density * NuTilde',integrand)


    end subroutine compute
    !******************************************************************************************









end module spalart_allmaras_laxfriedrichs
