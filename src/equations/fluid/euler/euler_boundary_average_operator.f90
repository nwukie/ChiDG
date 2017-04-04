module euler_boundary_average_operator
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, HALF
    use mod_fluid,              only: omega
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
        
        !
        ! Set operator name
        !
        call self%set_name("Euler Boundary Average Flux")

        !
        ! Set operator type
        !
        call self%set_operator_type("Boundary Advective Flux")

        !
        ! Set operator equations
        !
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
            u_a_m, v_a_m, w_a_m,                    &
            u_a_p, v_a_p, w_a_p,                    &
            energy_m,   energy_p,                   &
            enthalpy_m, enthalpy_p,                 &
            p_m,        p_p,                        &
            flux_1_m,   flux_2_m,   flux_3_m,       &
            flux_1_p,   flux_2_p,   flux_3_p,       &
            flux_1,     flux_2,     flux_3,         &
            integrand

        real(rk), allocatable, dimension(:) ::      &
            norm_1, norm_2, norm_3, r


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
        r = worker%coordinate('1','boundary') 
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / r
            mom2_p = mom2_p / r
        else if (worker%coordinate_system() == 'Cartesian') then

        else
            call chidg_signal(FATAL,"inlet, bad coordinate system")
        end if


        !
        ! Get normal vectors
        !
        norm_1 = worker%normal(1)
        norm_2 = worker%normal(2)
        norm_3 = worker%normal(3)


        
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
        ! Compute pressure and total enthalpy
        !
        p_m = worker%get_model_field_face('Pressure', 'value', 'face interior')
        p_p = worker%get_model_field_face('Pressure', 'value', 'face exterior')

        enthalpy_m = (energy_m + p_m)/density_m
        enthalpy_p = (energy_p + p_p)/density_p



        !=================================================
        ! mass flux
        !=================================================
        flux_1_m = (density_m * u_a_m)
        flux_2_m = (density_m * v_a_m)
        flux_3_m = (density_m * w_a_m)

        flux_1_p = (density_p * u_a_p)
        flux_2_p = (density_p * v_a_p)
        flux_3_p = (density_p * w_a_p)

        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)

        ! dot with normal vector
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)

        call worker%integrate_boundary('Density',integrand)


        !=================================================
        ! momentum-1 flux
        !=================================================
        flux_1_m = (mom1_m * u_a_m) + p_m
        flux_2_m = (mom1_m * v_a_m)
        flux_3_m = (mom1_m * w_a_m)

        flux_1_p = (mom1_p * u_a_p) + p_p
        flux_2_p = (mom1_p * v_a_p)
        flux_3_p = (mom1_p * w_a_p)

        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)


        ! dot with normal vector
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)

        call worker%integrate_boundary('Momentum-1',integrand)


        !=================================================
        ! momentum-2 flux
        !=================================================
        flux_1_m = (mom2_m * u_a_m)
        flux_2_m = (mom2_m * v_a_m) + p_m
        flux_3_m = (mom2_m * w_a_m)

        flux_1_p = (mom2_p * u_a_p)
        flux_2_p = (mom2_p * v_a_p) + p_p
        flux_3_p = (mom2_p * w_a_p)

        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)


        ! dot with normal vector
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)

        !
        ! Convert to tangential to angular momentum flux
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            integrand = integrand * r
        else if (worker%coordinate_system() == 'Cartesian') then

        else
            call chidg_signal(FATAL,"inlet, bad coordinate system")
        end if

        call worker%integrate_boundary('Momentum-2',integrand)


        !=================================================
        ! momentum-3 flux
        !=================================================
        flux_1_m = (mom3_m * u_a_m)
        flux_2_m = (mom3_m * v_a_m)
        flux_3_m = (mom3_m * w_a_m) + p_m

        flux_1_p = (mom3_p * u_a_p)
        flux_2_p = (mom3_p * v_a_p)
        flux_3_p = (mom3_p * w_a_p) + p_p


        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)


        ! dot with normal vector
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)

        call worker%integrate_boundary('Momentum-3',integrand)


        !=================================================
        ! energy flux
        !=================================================
        flux_1_m = (density_m * enthalpy_m * u_a_m)
        flux_2_m = (density_m * enthalpy_m * v_a_m)  +  omega*r*p_m
        flux_3_m = (density_m * enthalpy_m * w_a_m)

        flux_1_p = (density_p * enthalpy_p * u_a_p)
        flux_2_p = (density_p * enthalpy_p * v_a_p)  +  omega*r*p_p
        flux_3_p = (density_p * enthalpy_p * w_a_p)
        
        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)


        ! dot with normal vector
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)

        call worker%integrate_boundary('Energy',integrand)


    end subroutine compute
    !*********************************************************************************************************












end module euler_boundary_average_operator
