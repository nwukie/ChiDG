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


        type(AD_D), dimension(:), allocatable   ::                  &
            density_m, mom1_m, mom2_m, mom3_m, density_nutilde_m,   &
            density_p, mom1_p, mom2_p, mom3_p, density_nutilde_p,   &
            invdensity_m, u_m, v_m, w_m, T_m, un_m, c_m, diss_m,    &
            invdensity_p, u_p, v_p, w_p, T_p, un_p, c_p, diss_p,    &
            flux_avg_1, flux_avg_2, flux_avg_3, diff,               &
            flux_1, flux_2, flux_3, integrand

        real(rk),   dimension(:), allocatable   ::  &
            norm_1,  norm_2,  norm_3,               &
            unorm_1, unorm_2, unorm_3, r



        !
        ! Interpolate solution to quadrature nodes
        !
        density_m         = worker%get_primary_field_face('Density',           'value', 'face interior')
        mom1_m            = worker%get_primary_field_face('Momentum-1',        'value', 'face interior')
        mom2_m            = worker%get_primary_field_face('Momentum-2',        'value', 'face interior')
        mom3_m            = worker%get_primary_field_face('Momentum-3',        'value', 'face interior')
        density_nutilde_m = worker%get_primary_field_face('Density * NuTilde', 'value', 'face interior')


        density_p         = worker%get_primary_field_face('Density',           'value', 'face exterior')
        mom1_p            = worker%get_primary_field_face('Momentum-1',        'value', 'face exterior')
        mom2_p            = worker%get_primary_field_face('Momentum-2',        'value', 'face exterior')
        mom3_p            = worker%get_primary_field_face('Momentum-3',        'value', 'face exterior')
        density_nutilde_p = worker%get_primary_field_face('Density * NuTilde', 'value', 'face exterior')



        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / worker%coordinate('1','boundary')
            mom2_p = mom2_p / worker%coordinate('1','boundary')
        end if



        !
        ! Compute velocities
        !
        invdensity_m = ONE/density_m
        invdensity_p = ONE/density_p

        u_m = mom1_m*invdensity_m
        v_m = mom2_m*invdensity_m
        w_m = mom3_m*invdensity_m

        u_p = mom1_p*invdensity_p
        v_p = mom2_p*invdensity_p
        w_p = mom3_p*invdensity_p

        
        !
        ! Compute transport velocities
        !
        r = worker%coordinate('1','boundary') 
        v_m = v_m - omega*r
        v_p = v_p - omega*r


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
        ! Compute average flux and field difference.
        ! 
        flux_avg_1 = HALF*(u_m*density_nutilde_m  +  u_p*density_nutilde_p)
        flux_avg_2 = HALF*(v_m*density_nutilde_m  +  v_p*density_nutilde_p)
        flux_avg_3 = HALF*(w_m*density_nutilde_m  +  w_p*density_nutilde_p)



        un_m = u_m*unorm_1 + v_m*unorm_2 + w_m*unorm_3
        un_p = u_p*unorm_1 + v_p*unorm_2 + w_p*unorm_3


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

!        integrand = flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3
!        integrand = flux_1*norm_1*unorm_1 + flux_2*norm_2*unorm_2 + flux_3*norm_3*unorm_3

        integrand = flux_avg_1*norm_1 + flux_avg_2*norm_2 + flux_avg_3*norm_3

        integrand = integrand + max(abs(diss_m),abs(diss_p))*HALF*diff*norm_1*unorm_1
        integrand = integrand + max(abs(diss_m),abs(diss_p))*HALF*diff*norm_2*unorm_2
        integrand = integrand + max(abs(diss_m),abs(diss_p))*HALF*diff*norm_3*unorm_3


        !
        ! Integrate flux
        !
        call worker%integrate_boundary('Density * NuTilde',integrand)


    end subroutine compute
    !******************************************************************************************









end module spalart_allmaras_laxfriedrichs
