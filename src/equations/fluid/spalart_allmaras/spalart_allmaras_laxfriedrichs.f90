module spalart_allmaras_laxfriedrichs
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF
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
        call self%set_name("Spalart-Allmaras LaxFriedrichs Operator")

        ! Set operator type
        call self%set_operator_type("Boundary Advective Operator")

        ! Set operator equations
        call self%add_primary_field("Density * NuTilde")

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


        type(AD_D), dimension(:), allocatable   ::  &
            rho_m, rhou_m, rhov_m, rhow_m, rho_nutilde_m, invrho_m, u_m, v_m, w_m, T_m, un_m, c_m, diss_m, &
            rho_p, rhou_p, rhov_p, rhow_p, rho_nutilde_p, invrho_p, u_p, v_p, w_p, T_p, un_p, c_p, diss_p, &
            flux_avg_x, flux_avg_y, flux_avg_z, diff, &
            flux_x, flux_y, flux_z, integrand

        real(rk),   dimension(:), allocatable   ::  &
            normx, normy, normz, unormx, unormy, unormz


        !
        ! Interpolate solution to quadrature nodes
        !
        rho_m         = worker%get_primary_field_face('Density',           'value', 'face interior')
        rhou_m        = worker%get_primary_field_face('X-Momentum',        'value', 'face interior')
        rhov_m        = worker%get_primary_field_face('Y-Momentum',        'value', 'face interior')
        rhow_m        = worker%get_primary_field_face('Z-Momentum',        'value', 'face interior')
        rho_nutilde_m = worker%get_primary_field_face('Density * NuTilde', 'value', 'face interior')


        rho_p         = worker%get_primary_field_face('Density',           'value', 'face exterior')
        rhou_p        = worker%get_primary_field_face('X-Momentum',        'value', 'face exterior')
        rhov_p        = worker%get_primary_field_face('Y-Momentum',        'value', 'face exterior')
        rhow_p        = worker%get_primary_field_face('Z-Momentum',        'value', 'face exterior')
        rho_nutilde_p = worker%get_primary_field_face('Density * NuTilde', 'value', 'face exterior')


        T_m = worker%get_model_field_face('Temperature','value','face interior')
        T_p = worker%get_model_field_face('Temperature','value','face exterior')




        !
        ! Compute velocities
        !
        invrho_m = ONE/rho_m
        u_m = rhou_m*invrho_m
        v_m = rhov_m*invrho_m
        w_m = rhow_m*invrho_m

        invrho_p = ONE/rho_p
        u_p = rhou_p*invrho_p
        v_p = rhov_p*invrho_p
        w_p = rhow_p*invrho_p


        !
        ! Get normal vector
        !
        normx  = worker%normal(1)
        normy  = worker%normal(2)
        normz  = worker%normal(3)
        unormx = worker%unit_normal(1)
        unormy = worker%unit_normal(2)
        unormz = worker%unit_normal(3)


        !
        ! Compute average flux and field difference.
        ! 
        flux_avg_x = HALF*(u_m*rho_nutilde_m  +  u_p*rho_nutilde_p)
        flux_avg_y = HALF*(v_m*rho_nutilde_m  +  v_p*rho_nutilde_p)
        flux_avg_z = HALF*(w_m*rho_nutilde_m  +  w_p*rho_nutilde_p)



        un_m = u_m*unormx + v_m*unormy + w_m*unormz
        un_p = u_p*unormx + v_p*unormy + w_p*unormz

        c_m = sqrt(1.4_rk * 287.15_rk * T_m)
        c_p = sqrt(1.4_rk * 287.15_rk * T_p)

        diss_m = abs(un_m) + c_m
        diss_p = abs(un_p) + c_p

        !
        ! Compute Lax-Friedrichs upwind flux
        !
        diff = (rho_nutilde_m - rho_nutilde_p)
        flux_x = flux_avg_x + max(abs(diss_m),abs(diss_p))*HALF*diff
        flux_y = flux_avg_y + max(abs(diss_m),abs(diss_p))*HALF*diff
        flux_z = flux_avg_z + max(abs(diss_m),abs(diss_p))*HALF*diff

!        integrand = flux_x*normx + flux_y*normy + flux_z*normz
!        integrand = flux_x*normx*unormx + flux_y*normy*unormy + flux_z*normz*unormz

        integrand = flux_avg_x*normx + flux_avg_y*normy + flux_avg_z*normz

        integrand = integrand + max(abs(diss_m),abs(diss_p))*HALF*diff*normx*unormx
        integrand = integrand + max(abs(diss_m),abs(diss_p))*HALF*diff*normy*unormy
        integrand = integrand + max(abs(diss_m),abs(diss_p))*HALF*diff*normz*unormz


        !
        ! Integrate flux
        !
        call worker%integrate_boundary('Density * NuTilde',integrand)


    end subroutine compute
    !******************************************************************************************









end module spalart_allmaras_laxfriedrichs
