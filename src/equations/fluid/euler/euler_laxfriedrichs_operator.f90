module euler_laxfriedrichs_operator
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF
    use mod_fluid,              only: omega
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none

    private


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!
    !----------------------------------------------------------------------------------
    type, extends(operator_t), public :: euler_laxfriedrichs_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type euler_laxfriedrichs_operator_t
    !**********************************************************************************










contains

    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_laxfriedrichs_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name("Euler LaxFriedrichs Flux")

        ! Set operator type
        call self%set_operator_type("Boundary Advective Flux")

        ! Set operator equations
        call self%add_primary_field('Density'   )
        call self%add_primary_field('Momentum-1')
        call self%add_primary_field('Momentum-2')
        call self%add_primary_field('Momentum-3')
        call self%add_primary_field('Energy'    )

    end subroutine init
    !********************************************************************************









    !>  Compute Lax-Friedrichs upwind flux
    !!
    !!  Dissipation = -alpha(u_a_m - u_a_p)
    !!
    !!  Alpha is the maximum wave speed
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/16/2016
    !!
    !!------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(euler_laxfriedrichs_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::    &
            density_m,      density_p,              &
            mom1_m,         mom1_p,                 &
            mom2_m,         mom2_p,                 &
            mom3_m,         mom3_p,                 &
            u_a_m,          u_a_p,                  &
            v_a_m,          v_a_p,                  &
            w_a_m,          w_a_p,                  &
            energy_m,       energy_p,               &
            p_m,            p_p,                    &
            un_m,           un_p,                   &
            a_m,            a_p,                    &
            wave_m,         wave_p,                 &
            upwind,         wave,                   &
            integrand

        real(rk), allocatable, dimension(:)    ::   &
            norm_1,  norm_2,  norm_3,               &
            unorm_1, unorm_2, unorm_3,              &
            r

        real(rk) :: gam_m, gam_p


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


        !
        ! Get normal vectors
        !
        norm_1  = worker%normal(1)
        norm_2  = worker%normal(2)
        norm_3  = worker%normal(3)

        unorm_1 = worker%unit_normal(1)
        unorm_2 = worker%unit_normal(2)
        unorm_3 = worker%unit_normal(3)



!        !
!        ! Compute 1/rho
!        !
!        invdensity_m = ONE/density_m
!        invdensity_p = ONE/density_p
!
!
!        !
!        ! Compute velocity components
!        !
!        u_a_m = mom1_m*invdensity_m
!        v_m = mom2_m*invdensity_m
!        w_a_m = mom3_m*invdensity_m
!
!        u_a_p = mom1_p*invdensity_p
!        v_p = mom2_p*invdensity_p
!        w_a_p = mom3_p*invdensity_p



        !
        ! Fluid advection velocity
        !
        !r = worker%coordinate('1','boundary')
        !v_a_m = v_m - omega*r
        !v_a_p = v_p - omega*r
        u_a_m = worker%get_model_field_face('Advection Velocity-1', 'value', 'face interior')
        v_a_m = worker%get_model_field_face('Advection Velocity-2', 'value', 'face interior')
        w_a_m = worker%get_model_field_face('Advection Velocity-3', 'value', 'face interior')

        u_a_p = worker%get_model_field_face('Advection Velocity-1', 'value', 'face exterior')
        v_a_p = worker%get_model_field_face('Advection Velocity-2', 'value', 'face exterior')
        w_a_p = worker%get_model_field_face('Advection Velocity-3', 'value', 'face exterior')




        !
        ! Compute pressure and gamma
        !
        p_m = worker%get_model_field_face('Pressure','value','face interior')
        p_p = worker%get_model_field_face('Pressure','value','face exterior')
        gam_m = 1.4_rk
        gam_p = 1.4_rk


        !
        ! Compute normal velocities: dot-product vector projection along unit-normal direction
        !
        un_m = unorm_1*(u_a_m) + unorm_2*(v_a_m) + unorm_3*(w_a_m)
        un_p = unorm_1*(u_a_p) + unorm_2*(v_a_p) + unorm_3*(w_a_p)

        
        !
        ! Compute speed of sound
        !
        a_m = sqrt(abs(gam_m * p_m / density_m))
        a_p = sqrt(abs(gam_p * p_p / density_p))


        !
        ! Compute wave speeds
        !
        wave_m = abs(un_m) + a_m
        wave_p = abs(un_p) + a_p
        wave   = max(wave_m,wave_p)



        !===================================================
        ! mass flux
        !===================================================
        upwind = -wave*(density_p - density_m)

        integrand = HALF*(upwind*norm_1*unorm_1  +  upwind*norm_2*unorm_2  +  upwind*norm_3*unorm_3)

        call worker%integrate_boundary('Density',integrand)


        !===================================================
        ! momentum-1 flux
        !===================================================
        upwind = -wave*(mom1_p - mom1_m)

        integrand = HALF*(upwind*norm_1*unorm_1  +  upwind*norm_2*unorm_2  +  upwind*norm_3*unorm_3)

        call worker%integrate_boundary('Momentum-1',integrand)


        !===================================================
        ! momentum-2 flux
        !===================================================
        upwind = -wave*(mom2_p - mom2_m)

        integrand = HALF*(upwind*norm_1*unorm_1  +  upwind*norm_2*unorm_2  +  upwind*norm_3*unorm_3)

        !
        ! Convert to tangential to angular momentum flux
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            integrand = integrand * worker%coordinate('1','boundary')
        end if

        call worker%integrate_boundary('Momentum-2',integrand)

        !===================================================
        ! momentum-3 flux
        !===================================================
        upwind = -wave*(mom3_p - mom3_m)

        integrand = HALF*(upwind*norm_1*unorm_1  +  upwind*norm_2*unorm_2  +  upwind*norm_3*unorm_3)

        call worker%integrate_boundary('Momentum-3',integrand)

        !===================================================
        ! energy flux
        !===================================================
        upwind = -wave*(energy_p - energy_m)

        integrand = HALF*(upwind*norm_1*unorm_1  +  upwind*norm_2*unorm_2  +  upwind*norm_3*unorm_3)

        call worker%integrate_boundary('Energy',integrand)


    end subroutine compute
    !*******************************************************************************************













end module euler_laxfriedrichs_operator
