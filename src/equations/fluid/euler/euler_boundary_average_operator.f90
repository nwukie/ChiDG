module euler_boundary_average_operator
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
    type, extends(operator_t), public :: euler_boundary_average_operator_t

    contains
        procedure   :: init
        procedure   :: compute
    end type euler_boundary_average_operator_t
    !********************************************************************************


contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_boundary_average_operator_t),   intent(inout) :: self
        
        ! Set operator name
        call self%set_name("Euler Boundary Average Flux")

        ! Set operator type
        call self%set_operator_type("Boundary Advective Flux")

        ! Set operator equations
        call self%add_primary_field("Density"   )
        call self%add_primary_field("Momentum-1")
        call self%add_primary_field("Momentum-2")
        call self%add_primary_field("Momentum-3")
        call self%add_primary_field("Energy"    )

    end subroutine init
    !********************************************************************************



    !>  Boundary Flux routine for Euler
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!-------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(euler_boundary_average_operator_t),   intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable,    dimension(:) :: &
            density_m,  density_p,                  &
            mom1_m,     mom1_p,                     &
            mom2_m,     mom2_p,                     &
            mom3_m,     mom3_p,                     &
            energy_m,   energy_p,                   &
            p_m,        p_p,                        &
            enthalpy_m, enthalpy_p,                 &
            u_m,        v_m,        w_m,            &
            u_p,        v_p,        w_p,            &
            flux_1_m,   flux_2_m,   flux_3_m,       &
            flux_1_p,   flux_2_p,   flux_3_p,       &
            invdensity_m,   invdensity_p, r


        ! Interpolate solution to quadrature nodes
        density_m = worker%get_field('Density'   , 'value', 'face interior')
        density_p = worker%get_field('Density'   , 'value', 'face exterior')

        mom1_m    = worker%get_field('Momentum-1', 'value', 'face interior')
        mom1_p    = worker%get_field('Momentum-1', 'value', 'face exterior')

        mom2_m    = worker%get_field('Momentum-2', 'value', 'face interior')
        mom2_p    = worker%get_field('Momentum-2', 'value', 'face exterior')

        mom3_m    = worker%get_field('Momentum-3', 'value', 'face interior')
        mom3_p    = worker%get_field('Momentum-3', 'value', 'face exterior')

        energy_m  = worker%get_field('Energy'    , 'value', 'face interior')
        energy_p  = worker%get_field('Energy'    , 'value', 'face exterior')


        ! Account for cylindrical. Get tangential momentum from angular momentum.
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','face interior') 
            mom2_m = mom2_m / r
            mom2_p = mom2_p / r
        end if


        ! Compute velocity
        invdensity_m = ONE/density_m
        invdensity_p = ONE/density_p

        u_m = mom1_m*invdensity_m
        v_m = mom2_m*invdensity_m
        w_m = mom3_m*invdensity_m

        u_p = mom1_p*invdensity_p
        v_p = mom2_p*invdensity_p
        w_p = mom3_p*invdensity_p


        ! Compute pressure and total enthalpy
        p_m = worker%get_field('Pressure', 'value', 'face interior')
        p_p = worker%get_field('Pressure', 'value', 'face exterior')

        enthalpy_m = (energy_m + p_m)*invdensity_m
        enthalpy_p = (energy_p + p_p)*invdensity_p


        !================================
        !           Mass flux
        !================================
        flux_1_m = (density_m * u_m)
        flux_2_m = (density_m * v_m)
        flux_3_m = (density_m * w_m)

        flux_1_p = (density_p * u_p)
        flux_2_p = (density_p * v_p)
        flux_3_p = (density_p * w_p)

        call worker%integrate_boundary_average('Density','Advection',       &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)


        !================================
        !           Momentum-1 
        !================================
        flux_1_m = (mom1_m * u_m) + p_m
        flux_2_m = (mom1_m * v_m)
        flux_3_m = (mom1_m * w_m)

        flux_1_p = (mom1_p * u_p) + p_p
        flux_2_p = (mom1_p * v_p)
        flux_3_p = (mom1_p * w_p)

        call worker%integrate_boundary_average('Momentum-1','Advection',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)

        !================================
        !           Momentum-2
        !================================
        flux_1_m = (mom2_m * u_m)
        flux_2_m = (mom2_m * v_m) + p_m
        flux_3_m = (mom2_m * w_m)
                            
        flux_1_p = (mom2_p * u_p)
        flux_2_p = (mom2_p * v_p) + p_p
        flux_3_p = (mom2_p * w_p)

        ! Convert to tangential to angular momentum flux
        if (worker%coordinate_system() == 'Cylindrical') then
            flux_1_m = flux_1_m * r
            flux_2_m = flux_2_m * r
            flux_3_m = flux_3_m * r

            flux_1_p = flux_1_p * r
            flux_2_p = flux_2_p * r
            flux_3_p = flux_3_p * r
        end if

        call worker%integrate_boundary_average('Momentum-2','Advection',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)


        !================================
        !           Momentum-3
        !================================
        flux_1_m = (mom3_m * u_m)
        flux_2_m = (mom3_m * v_m)
        flux_3_m = (mom3_m * w_m) + p_m
                    
        flux_1_p = (mom3_p * u_p)
        flux_2_p = (mom3_p * v_p)
        flux_3_p = (mom3_p * w_p) + p_p

        call worker%integrate_boundary_average('Momentum-3','Advection',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)


        !================================
        !           Energy
        !================================
        flux_1_m = (density_m * enthalpy_m * u_m)
        flux_2_m = (density_m * enthalpy_m * v_m)
        flux_3_m = (density_m * enthalpy_m * w_m)

        flux_1_p = (density_p * enthalpy_p * u_p)
        flux_2_p = (density_p * enthalpy_p * v_p)
        flux_3_p = (density_p * enthalpy_p * w_p)

        call worker%integrate_boundary_average('Energy','Advection',        &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)


    end subroutine compute
    !**********************************************************************************












end module euler_boundary_average_operator
