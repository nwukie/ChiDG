module fluid_viscous_boundary_average_operator
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, HALF

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
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
    type, extends(operator_t), public :: fluid_viscous_boundary_average_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type fluid_viscous_boundary_average_operator_t
    !********************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(fluid_viscous_boundary_average_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name('Fluid Viscous Boundary Average Operator')

        !
        ! Set operator type
        !
        call self%set_operator_type('Boundary Diffusive Flux')

        !
        ! Set operator equations
        !
        call self%add_primary_field('Density'   )
        call self%add_primary_field('Momentum-1')
        call self%add_primary_field('Momentum-2')
        call self%add_primary_field('Momentum-3')
        call self%add_primary_field('Energy'    )

    end subroutine init
    !********************************************************************************



    !>  Boundary Flux routine for Euler
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!-------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(fluid_viscous_boundary_average_operator_t),   intent(inout)   :: self
        type(chidg_worker_t),                               intent(inout)   :: worker
        class(properties_t),                                intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::                        &
            density_m, mom1_m, mom2_m, mom3_m, energy_m,                &
            density_p, mom1_p, mom2_p, mom3_p, energy_p,                &
            u_m, v_m, w_m, invdensity_m,                                &
            u_p, v_p, w_p, invdensity_p,                                &
            k_m, k_p, k_l_m, k_l_p, k_t_m, k_t_p,                       &
            grad1_T_m, grad2_T_m, grad3_T_m,                            &
            grad1_T_p, grad2_T_p, grad3_T_p,                            &
            tau_11_m, tau_22_m, tau_33_m, tau_12_m, tau_13_m, tau_23_m, &
            tau_11_p, tau_22_p, tau_33_p, tau_12_p, tau_13_p, tau_23_p, &
            flux_1_m, flux_2_m, flux_3_m,                               &
            flux_1_p, flux_2_p, flux_3_p,                               &
            flux_1, flux_2, flux_3, integrand


        real(rk), allocatable, dimension(:) ::      &
            norm_1, norm_2, norm_3


        !
        ! Interpolate solution to quadrature nodes
        !
        density_m = worker%get_primary_field_face('Density'   , 'value', 'face interior')
        density_p = worker%get_primary_field_face('Density'   , 'value', 'face exterior')

        mom1_m    = worker%get_primary_field_face('Momentum-1', 'value', 'face interior')
        mom1_p    = worker%get_primary_field_face('Momentum-1', 'value', 'face exterior')

        mom2_m    = worker%get_primary_field_face('Momentum-2', 'value', 'face interior')
        mom2_p    = worker%get_primary_field_face('Momentum-2', 'value', 'face exterior')

        mom3_m    = worker%get_primary_field_face('Momentum-3', 'value', 'face interior')
        mom3_p    = worker%get_primary_field_face('Momentum-3', 'value', 'face exterior')

        energy_m  = worker%get_primary_field_face('Energy'    , 'value', 'face interior')
        energy_p  = worker%get_primary_field_face('Energy'    , 'value', 'face exterior')

        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / worker%coordinate('1','boundary')
            mom2_p = mom2_p / worker%coordinate('1','boundary')
        end if


        invdensity_m = ONE/density_m
        invdensity_p = ONE/density_p


        !
        ! Get normal vector
        !
        norm_1 = worker%normal(1)
        norm_2 = worker%normal(2)
        norm_3 = worker%normal(3)




        !
        ! Get model fields:
        !   Pressure
        !   Temperature
        !   Viscosity
        !   Second Coefficient of Viscosity
        !   Thermal Conductivity
        !
        k_l_m = worker%get_model_field_face('Laminar Thermal Conductivity',   'value', 'face interior')
        k_l_p = worker%get_model_field_face('Laminar Thermal Conductivity',   'value', 'face exterior')

        k_t_m = worker%get_model_field_face('Turbulent Thermal Conductivity', 'value', 'face interior')
        k_t_p = worker%get_model_field_face('Turbulent Thermal Conductivity', 'value', 'face exterior')



        !
        ! Compute effective viscosities, conductivity. Laminar + Turbulent
        !
        k_m = k_l_m + k_t_m
        k_p = k_l_p + k_t_p



        !
        ! Compute velocities
        !
        u_m = mom1_m/density_m
        v_m = mom2_m/density_m
        w_m = mom3_m/density_m

        u_p = mom1_p/density_p
        v_p = mom2_p/density_p
        w_p = mom3_p/density_p



        !
        ! get temperature gradient
        !
        grad1_T_m = worker%get_model_field_face('Temperature Gradient - 1', 'value', 'face interior')
        grad2_T_m = worker%get_model_field_face('Temperature Gradient - 2', 'value', 'face interior')
        grad3_T_m = worker%get_model_field_face('Temperature Gradient - 3', 'value', 'face interior')

        grad1_T_p = worker%get_model_field_face('Temperature Gradient - 1', 'value', 'face exterior')
        grad2_T_p = worker%get_model_field_face('Temperature Gradient - 2', 'value', 'face exterior')
        grad3_T_p = worker%get_model_field_face('Temperature Gradient - 3', 'value', 'face exterior')




        !
        ! get shear stress components
        !
        tau_11_m = worker%get_model_field_face('Shear-11', 'value', 'face interior')
        tau_22_m = worker%get_model_field_face('Shear-22', 'value', 'face interior')
        tau_33_m = worker%get_model_field_face('Shear-33', 'value', 'face interior')

        tau_12_m = worker%get_model_field_face('Shear-12', 'value', 'face interior')
        tau_13_m = worker%get_model_field_face('Shear-13', 'value', 'face interior')
        tau_23_m = worker%get_model_field_face('Shear-23', 'value', 'face interior')


        tau_11_p = worker%get_model_field_face('Shear-11', 'value', 'face exterior')
        tau_22_p = worker%get_model_field_face('Shear-22', 'value', 'face exterior')
        tau_33_p = worker%get_model_field_face('Shear-33', 'value', 'face exterior')

        tau_12_p = worker%get_model_field_face('Shear-12', 'value', 'face exterior')
        tau_13_p = worker%get_model_field_face('Shear-13', 'value', 'face exterior')
        tau_23_p = worker%get_model_field_face('Shear-23', 'value', 'face exterior')



        !----------------------------------
        !            mass flux
        !----------------------------------


        !----------------------------------
        !         momentum-1 flux
        !----------------------------------
        flux_1_m = -tau_11_m
        flux_2_m = -tau_12_m
        flux_3_m = -tau_13_m

        flux_1_p = -tau_11_p
        flux_2_p = -tau_12_p
        flux_3_p = -tau_13_p

        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)


        ! dot with normal vector
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)

        call worker%integrate_boundary('Momentum-1',integrand)


        !----------------------------------
        !         momentum-2 flux
        !----------------------------------
        flux_1_m = -tau_12_m
        flux_2_m = -tau_22_m
        flux_3_m = -tau_23_m

        flux_1_p = -tau_12_p
        flux_2_p = -tau_22_p
        flux_3_p = -tau_23_p

        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)


        ! dot with normal vector
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)

        !
        ! Convert to tangential to angular momentum flux
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            integrand = integrand * worker%coordinate('1','boundary')
        end if


        call worker%integrate_boundary('Momentum-2',integrand)


        !----------------------------------
        !         momentum-3 flux
        !----------------------------------
        flux_1_m = -tau_13_m
        flux_2_m = -tau_23_m
        flux_3_m = -tau_33_m

        flux_1_p = -tau_13_p
        flux_2_p = -tau_23_p
        flux_3_p = -tau_33_p


        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)


        ! dot with normal vector
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)

        call worker%integrate_boundary('Momentum-3',integrand)


        !----------------------------------
        !           energy flux
        !----------------------------------
        flux_1_m = -k_m*grad1_T_m  -  (u_m*tau_11_m + v_m*tau_12_m + w_m*tau_13_m)
        flux_2_m = -k_m*grad2_T_m  -  (u_m*tau_12_m + v_m*tau_22_m + w_m*tau_23_m)
        flux_3_m = -k_m*grad3_T_m  -  (u_m*tau_13_m + v_m*tau_23_m + w_m*tau_33_m)

        flux_1_p = -k_p*grad1_T_p  -  (u_p*tau_11_p + v_p*tau_12_p + w_p*tau_13_p)
        flux_2_p = -k_p*grad2_T_p  -  (u_p*tau_12_p + v_p*tau_22_p + w_p*tau_23_p)
        flux_3_p = -k_p*grad3_T_p  -  (u_p*tau_13_p + v_p*tau_23_p + w_p*tau_33_p)


        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)


        ! dot with normal vector
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)

        call worker%integrate_boundary('Energy',integrand)


    end subroutine compute
    !*********************************************************************************************************












end module fluid_viscous_boundary_average_operator
