module euler_ale_laxfriedrichs_operator
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF
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
    type, extends(operator_t), public :: euler_ale_laxfriedrichs_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type euler_ale_laxfriedrichs_operator_t
    !**********************************************************************************










contains

    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_ale_laxfriedrichs_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name("Euler ALE LaxFriedrichs Flux")

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
        class(euler_ale_laxfriedrichs_operator_t),  intent(inout)   :: self
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

        real(rk), allocatable, dimension(:) ::      &
           u_grid, v_grid, w_grid, det_jacobian_grid


        real(rk), allocatable, dimension(:,:,:) ::      &
            jacobian_grid


!        print *, 'LaxFriedrichs Flux'
        u_grid = worker%get_grid_velocity_face("u_grid")
        v_grid = worker%get_grid_velocity_face("v_grid")
        w_grid = worker%get_grid_velocity_face("w_grid")

        jacobian_grid = worker%get_inv_jacobian_grid_face()
        det_jacobian_grid = worker%get_det_jacobian_grid_face()




        irho  = prop%get_primary_field_index("Density"   )
        irhou = prop%get_primary_field_index("Momentum-1")
        irhov = prop%get_primary_field_index("Momentum-2")
        irhow = prop%get_primary_field_index("Momentum-3")
        irhoE = prop%get_primary_field_index("Energy"    )



        !
        ! Interpolate solution to quadrature nodes
        !
        rho_m  = worker%get_primary_field_face('Density'   , 'value', 'face interior')
        rho_p  = worker%get_primary_field_face('Density'   , 'value', 'face exterior')

        rhou_m = worker%get_primary_field_face('Momentum-1', 'value', 'face interior')
        rhou_p = worker%get_primary_field_face('Momentum-1', 'value', 'face exterior')

        rhov_m = worker%get_primary_field_face('Momentum-2', 'value', 'face interior')
        rhov_p = worker%get_primary_field_face('Momentum-2', 'value', 'face exterior')

        rhow_m = worker%get_primary_field_face('Momentum-3', 'value', 'face interior')
        rhow_p = worker%get_primary_field_face('Momentum-3', 'value', 'face exterior')

        rhoE_m = worker%get_primary_field_face('Energy'    , 'value', 'face interior')
        rhoE_p = worker%get_primary_field_face('Energy'    , 'value', 'face exterior')



        rho_m = rho_m/det_jacobian_grid
        rhou_m = rhou_m/det_jacobian_grid
        rhov_m = rhov_m/det_jacobian_grid
        rhow_m = rhow_m/det_jacobian_grid
        rhoE_m = rhoE_m/det_jacobian_grid

        rho_p = rho_p/det_jacobian_grid
        rhou_p = rhou_p/det_jacobian_grid
        rhov_p = rhov_p/det_jacobian_grid
        rhow_p = rhow_p/det_jacobian_grid
        rhoE_p = rhoE_p/det_jacobian_grid




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
        !p_m = worker%get_model_field_face('Pressure','value','face interior')
        !p_p = worker%get_model_field_face('Pressure','value','face exterior')
        p_m = (1.4_rk-ONE)*(rhoE_m - HALF*(rhou_m**TWO+rhov_m**TWO+rhow_m**TWO)/rho_m)
        p_p = (1.4_rk-ONE)*(rhoE_p - HALF*(rhou_p**TWO+rhov_p**TWO+rhow_p**TWO)/rho_p)


        gam_m = 1.4_rk
        gam_p = 1.4_rk


        !
        ! Compute normal velocities: dot-product vector projection along unit-normal direction
        !
        un_m = unormx*(rhou_m/rho_m-u_grid) + unormy*(rhov_m/rho_m-v_grid) + unormz*(rhow_m/rho_m-w_grid)
        un_p = unormx*(rhou_p/rho_p-u_grid) + unormy*(rhov_p/rho_p-v_grid) + unormz*(rhow_p/rho_p-w_grid)

        
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


!        print *, 'normx'
!        print *, normx
!        print *, 'normy'
!        print *, normy
!        print *, 'normz'
!        print *, normz



        allocate(norm_mag(size(normx)))
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

        call worker%integrate_boundary('Momentum-1',integrand)


        !================================
        !       Y-MOMENTUM FLUX
        !================================
        upwind = -wave*(rhov_p - rhov_m)

        integrand = HALF * upwind * norm_mag

        call worker%integrate_boundary('Momentum-2',integrand)

        !================================
        !       Z-MOMENTUM FLUX
        !================================
        upwind = -wave*(rhow_p - rhow_m)

        integrand = HALF * upwind * norm_mag

        call worker%integrate_boundary('Momentum-3',integrand)

        !================================
        !          ENERGY FLUX
        !================================
        upwind = -wave*(rhoE_p - rhoE_m)

        integrand = HALF * upwind * norm_mag

        call worker%integrate_boundary('Energy',integrand)


    end subroutine compute
    !*******************************************************************************************













end module euler_ale_laxfriedrichs_operator
