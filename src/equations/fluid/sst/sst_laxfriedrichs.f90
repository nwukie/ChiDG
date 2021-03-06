module sst_laxfriedrichs
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,THREE, HALF
    use mod_fluid,              only: omega, gam, Rgas
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none


    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: sst_laxfriedrichs_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type sst_laxfriedrichs_operator_t
    !*****************************************************************************************

contains



    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(sst_laxfriedrichs_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('SST LaxFriedrichs Operator')

        ! Set operator type
        call self%set_operator_type('Boundary Advective Operator')

        ! Set operator equations
        call self%add_primary_field('Density * Omega')
        call self%add_primary_field('Density * k')
    end subroutine init
    !*****************************************************************************************



    !>
    !! 
    !!
    !! @author  Eric M. Wolf
    !! @date    02/01/2018 
    !!
    !--------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(sst_laxfriedrichs_operator_t), intent(inout)   :: self
        type(chidg_worker_t),                             intent(inout)   :: worker
        class(properties_t),                              intent(inout)   :: prop


        type(AD_D), dimension(:), allocatable   ::      &
            density_m, density_p,                       &
            mom1_m, mom2_m, mom3_m,                     &
            mom1_p, mom2_p, mom3_p,                     &
            density_omega_m, density_omega_p,       &
            density_k_m, density_k_p,       &
            k_m,        k_p,                            &
            invdensity_m, invdensity_p,                 &
            u_m, v_m, w_m, T_m, un_m, c_m, wavespeed_m, &
            u_p, v_p, w_p, T_p, un_p, c_p, wavespeed_p, &
            dissipation, unorm_1, unorm_2, unorm_3, grid_vel_n, r

        real(rk),   dimension(:,:), allocatable :: grid_vel



        !
        ! Interpolate solution to quadrature nodes
        !
        density_m         = worker%get_field('Density',           'value', 'face interior')
        density_p         = worker%get_field('Density',           'value', 'face exterior')

        mom1_m            = worker%get_field('Momentum-1',        'value', 'face interior')
        mom1_p            = worker%get_field('Momentum-1',        'value', 'face exterior')

        mom2_m            = worker%get_field('Momentum-2',        'value', 'face interior')
        mom2_p            = worker%get_field('Momentum-2',        'value', 'face exterior')

        mom3_m            = worker%get_field('Momentum-3',        'value', 'face interior')
        mom3_p            = worker%get_field('Momentum-3',        'value', 'face exterior')

        density_omega_m = worker%get_field('Density * Omega', 'value', 'face interior')
        density_omega_p = worker%get_field('Density * Omega', 'value', 'face exterior')

        density_k_m = worker%get_field('Density * k', 'value', 'face interior')
        density_k_p = worker%get_field('Density * k', 'value', 'face exterior')

        k_m = worker%get_field('k', 'value', 'face interior')
        k_p = worker%get_field('k', 'value', 'face exterior')
                !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','face interior') 
            mom2_m = mom2_m / r
            mom2_p = mom2_p / r
        end if

        
        !
        ! Get fluid advection velocity
        !
        invdensity_m = ONE/density_m
        u_m = mom1_m*invdensity_m
        v_m = mom2_m*invdensity_m
        w_m = mom3_m*invdensity_m

        invdensity_p = ONE/density_p
        u_p = mom1_p*invdensity_p
        v_p = mom2_p*invdensity_p
        w_p = mom3_p*invdensity_p
        


        !
        ! Get normal vector
        !
        unorm_1 = worker%unit_normal_ale(1)
        unorm_2 = worker%unit_normal_ale(2)
        unorm_3 = worker%unit_normal_ale(3)



        !
        ! Compute maximum wave speed
        !
        grid_vel = worker%get_grid_velocity_face('face interior')
        un_m = u_m*unorm_1 + v_m*unorm_2 + w_m*unorm_3
        un_p = u_p*unorm_1 + v_p*unorm_2 + w_p*unorm_3
        grid_vel_n = grid_vel(:,1)*unorm_1  +  grid_vel(:,2)*unorm_2  +  grid_vel(:,3)*unorm_3


        T_m = worker%get_field('Temperature','value','face interior')
        T_p = worker%get_field('Temperature','value','face exterior')
        c_m = sqrt(gam * Rgas * T_m + (TWO/THREE)*1.4_rk*k_m)
        c_p = sqrt(gam * Rgas * T_p + (TWO/THREE)*1.4_rk*k_p)

        wavespeed_m = abs(un_m) + c_m  + abs(grid_vel_n)
        wavespeed_p = abs(un_p) + c_p  + abs(grid_vel_n)

        wavespeed_m = wavespeed_m
        wavespeed_p = wavespeed_p

        !
        ! Compute Lax-Friedrichs upwind flux
        !
        dissipation = HALF*max(abs(wavespeed_m),abs(wavespeed_p))*(density_omega_m - density_omega_p)


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_upwind('Density * Omega',dissipation)

        !
        ! Compute Lax-Friedrichs upwind flux
        !
        dissipation = HALF*max(abs(wavespeed_m),abs(wavespeed_p))*(density_k_m - density_k_p)


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_upwind('Density * k',dissipation)

        
    end subroutine compute
    !******************************************************************************************









end module sst_laxfriedrichs
