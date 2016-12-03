module euler_laxfriedrichs_operator
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO,HALF
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
        call self%add_primary_field("Density"   )
        call self%add_primary_field("X-Momentum")
        call self%add_primary_field("Y-Momentum")
        call self%add_primary_field("Z-Momentum")
        call self%add_primary_field("Energy"    )

    end subroutine init
    !********************************************************************************









    !>  Compute Lax-Friedrichs upwind flux
    !!
    !!  Dissipation = -alpha(u_m - u_p)
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

        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) :: &
            rho_m,      rho_p,                   &
            rhou_m,     rhou_p,                  &
            rhov_m,     rhov_p,                  &
            rhow_m,     rhow_p,                  &
            rhoe_m,     rhoe_p,                  &
            p_m,        p_p,                     &
            un_m,       un_p,                    &
            a_m,        a_p,                     &
            wave_m,     wave_p,                  &
            upwind,     wave,                    &
            integrand

        real(rk), allocatable, dimension(:)    :: &
            norm_mag, normx, normy, normz, unormx, unormy, unormz

        real(rk) :: gam_m, gam_p


        irho  = prop%get_primary_field_index("Density"   )
        irhou = prop%get_primary_field_index("X-Momentum")
        irhov = prop%get_primary_field_index("Y-Momentum")
        irhow = prop%get_primary_field_index("Z-Momentum")
        irhoE = prop%get_primary_field_index("Energy"    )



        !
        ! Interpolate solution to quadrature nodes
        !
        rho_m  = worker%get_primary_field_face('Density'   , 'value', 'face interior')
        rho_p  = worker%get_primary_field_face('Density'   , 'value', 'face exterior')

        rhou_m = worker%get_primary_field_face('X-Momentum', 'value', 'face interior')
        rhou_p = worker%get_primary_field_face('X-Momentum', 'value', 'face exterior')

        rhov_m = worker%get_primary_field_face('Y-Momentum', 'value', 'face interior')
        rhov_p = worker%get_primary_field_face('Y-Momentum', 'value', 'face exterior')

        rhow_m = worker%get_primary_field_face('Z-Momentum', 'value', 'face interior')
        rhow_p = worker%get_primary_field_face('Z-Momentum', 'value', 'face exterior')

        rhoE_m = worker%get_primary_field_face('Energy'    , 'value', 'face interior')
        rhoE_p = worker%get_primary_field_face('Energy'    , 'value', 'face exterior')






        normx  = worker%normal(1)
        normy  = worker%normal(2)
        normz  = worker%normal(3)

        unormx = worker%unit_normal(1)
        unormy = worker%unit_normal(2)
        unormz = worker%unit_normal(3)




        !
        ! Compute pressure and gamma
        !
        !p_m = prop%fluid%compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m)
        !p_p = prop%fluid%compute_pressure(rho_p,rhou_p,rhov_p,rhow_p,rhoE_p)
        !gam_m = prop%fluid%compute_gamma(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m)
        !gam_p = prop%fluid%compute_gamma(rho_p,rhou_p,rhov_p,rhow_p,rhoE_p)
        p_m = worker%get_model_field_face('Pressure','value','face interior')
        p_p = worker%get_model_field_face('Pressure','value','face exterior')
        gam_m = 1.4_rk
        gam_p = 1.4_rk


        !
        ! Compute normal velocities: dot-product vector projection along unit-normal direction
        !
        un_m = unormx*(rhou_m/rho_m) + unormy*(rhov_m/rho_m) + unormz*(rhow_m/rho_m)
        un_p = unormx*(rhou_p/rho_p) + unormy*(rhov_p/rho_p) + unormz*(rhow_p/rho_p)

        
        !
        ! Compute speed of sound
        !
        a_m = sqrt(abs(gam_m * p_m / rho_m))
        a_p = sqrt(abs(gam_p * p_p / rho_p))


        !
        ! Compute wave speeds
        !
        wave_m = abs(un_m) + a_m
        wave_p = abs(un_p) + a_p
        wave   = max(wave_m,wave_p)


        norm_mag = sqrt(normx**TWO + normy**TWO + normz**TWO)

        !================================
        !       MASS FLUX
        !================================
        upwind = -wave*(rho_p - rho_m)

        integrand = HALF * upwind * norm_mag

        call worker%integrate_boundary('Density',integrand)


        !================================
        !       X-MOMENTUM FLUX
        !================================
        upwind = -wave*(rhou_p - rhou_m)

        integrand = HALF * upwind * norm_mag

        call worker%integrate_boundary('X-Momentum',integrand)


        !================================
        !       Y-MOMENTUM FLUX
        !================================
        upwind = -wave*(rhov_p - rhov_m)

        integrand = HALF * upwind * norm_mag

        call worker%integrate_boundary('Y-Momentum',integrand)

        !================================
        !       Z-MOMENTUM FLUX
        !================================
        upwind = -wave*(rhow_p - rhow_m)

        integrand = HALF * upwind * norm_mag

        call worker%integrate_boundary('Z-Momentum',integrand)

        !================================
        !          ENERGY FLUX
        !================================
        upwind = -wave*(rhoE_p - rhoE_m)

        integrand = HALF * upwind * norm_mag

        call worker%integrate_boundary('Energy',integrand)


    end subroutine compute
    !*******************************************************************************************













end module euler_laxfriedrichs_operator
