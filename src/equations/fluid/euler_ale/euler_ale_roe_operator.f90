module euler_ale_roe_operator
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF
    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    use ieee_arithmetic,        only: ieee_is_nan
    implicit none


    
    !> Implementation of Roe's approximate Riemann solver.
    !!
    !! The formulation used here is from the reference:
    !!   J. Blazek,"Computational Fluid Dynamics: Principles and Applications"
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: euler_ale_roe_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type euler_ale_roe_operator_t
    !*******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_ale_roe_operator_t),   intent(inout)    :: self

        ! Set operator name
        call self%set_name("Euler ALE Roe Flux")

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




    !> Compute Roe approximate Riemann upwind flux
    !! 
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!---------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(euler_ale_roe_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                      &
            rho_m,      rho_p,                                          &
            rhou_m,     rhou_p,                                         &
            rhov_m,     rhov_p,                                         &
            rhow_m,     rhow_p,                                         &
            rhoe_m,     rhoe_p,                                         &
            rhor_p,     rhot_p,                                         &
            p_m,        p_p,                                            &
            un_m,       un_p,                                           &
            a_m,        a_p,                                            &
            H_m,        H_p,                                            &
            rtil, util, vtil, wtil, vmagtil, Htil, ctil, qtil2,         &
            integrand,  upwind,     wave,                               &
            C1,  C2_a, C2_b,  C3,                                       &
            u_m, v_m, w_m,                                              &
            u_p, v_p, w_p,                                              &
            vmag_p, vmag_m,                                             &
            delr,   delp,   delvmag, delu, delv, delw,                  &
            lamda1, lamda2, lamda3,                                     &
            sqrt_rhom, sqrt_rhop, sqrt_rhom_plus_rhop, ctil2, invrho_m, invrho_p

        real(rk), allocatable, dimension(:) :: &
            normx, normy, normz, unormx, unormy, unormz

        real(rk), allocatable, dimension(:) ::      &
           u_grid, v_grid, w_grid, det_jacobian_grid


        real(rk), allocatable, dimension(:,:,:) ::      &
            jacobian_grid


        real(rk) :: eps, gam_m, gam_p

        u_grid = worker%get_grid_velocity_face("u_grid",'face interior')
        v_grid = worker%get_grid_velocity_face("v_grid",'face interior')
        w_grid = worker%get_grid_velocity_face("w_grid",'face interior')

        jacobian_grid = worker%get_inv_jacobian_grid_face('face interior')
        det_jacobian_grid = worker%get_det_jacobian_grid_face('value', 'face interior')



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
        !p_m = worker%get_model_field_face('Pressure', 'value', 'face interior')
        !p_p = worker%get_model_field_face('Pressure', 'value', 'face exterior')
        gam_m = 1.4_rk
        gam_p = 1.4_rk

        p_m = (1.4_rk-ONE)*(rhoE_m - HALF*(rhou_m**TWO+rhov_m**TWO+rhow_m**TWO)/rho_m)
        p_p = (1.4_rk-ONE)*(rhoE_p - HALF*(rhou_p**TWO+rhov_p**TWO+rhow_p**TWO)/rho_p)

        invrho_m = ONE/rho_m
        invrho_p = ONE/rho_p

        !
        ! Compute enthalpy
        !
        H_m = (rhoE_m + p_m)*invrho_m
        H_p = (rhoE_p + p_p)*invrho_p


        !
        ! Compute velocity components
        !
        u_m = rhou_m*invrho_m 
        v_m = rhov_m*invrho_m 
        w_m = rhow_m*invrho_m 
        vmag_m = u_m*unormx + v_m*unormy + w_m*unormz

        u_p = rhou_p*invrho_p 
        v_p = rhov_p*invrho_p 
        w_p = rhow_p*invrho_p 
        vmag_p = u_p*unormx + v_p*unormy + w_p*unormz


        !
        ! Compute Roe-averaged variables
        !
        sqrt_rhom = sqrt(rho_m)
        sqrt_rhop = sqrt(rho_p)
        sqrt_rhom_plus_rhop = sqrt_rhom + sqrt_rhop
        rtil =  sqrt(rho_p * rho_m)                                       ! Roe-averaged density
        util = (sqrt_rhom*u_m + sqrt_rhop*u_p) / (sqrt_rhom_plus_rhop)    ! Roe-averaged u-velocity
        vtil = (sqrt_rhom*v_m + sqrt_rhop*v_p) / (sqrt_rhom_plus_rhop)    ! Roe-averaged v-velocity
        wtil = (sqrt_rhom*w_m + sqrt_rhop*w_p) / (sqrt_rhom_plus_rhop)    ! Roe-averaged w-velocity
        Htil = (sqrt_rhom*H_m + sqrt_rhop*H_p) / (sqrt_rhom_plus_rhop)    ! Roe-averaged Enthalpy

        vmagtil = util*unormx + vtil*unormy + wtil*unormz  ! Magnitude of Roe-averaged velocity in the face normal direction
        qtil2   = util**TWO + vtil**TWO + wtil**TWO


        !& HARDCODED GAMMA
        ctil = sqrt((1.4_rk - ONE)*(Htil - HALF*qtil2))                   ! Roe-averaged speed of sound
        ctil2 = ctil**TWO



        !
        ! Compute jump terms
        !
        delr    = (rho_m - rho_p)
        delu    = (u_m - u_p)
        delv    = (v_m - v_p)
        delw    = (w_m - w_p)
        delvmag = (vmag_m - vmag_p)
        delp    = (p_m - p_p)


        !
        ! Limit wave speeds for entropy fix
        !
        lamda1 = abs(vmagtil - ctil) + sqrt(u_grid**TWO+v_grid**TWO+w_grid**TWO)
        lamda2 = abs(vmagtil) + sqrt(u_grid**TWO+v_grid**TWO+w_grid**TWO)
        lamda3 = abs(vmagtil + ctil) + sqrt(u_grid**TWO+v_grid**TWO+w_grid**TWO)

        eps = 0.01_rk
        where ( (-eps*ctil < lamda1) .and. (lamda1 < eps*ctil) )
            lamda1 = HALF*(eps*ctil + lamda1*lamda1/(eps*ctil))
        else where
            lamda1 = sqrt(lamda1*lamda1)
        end where

        where ( (-eps*ctil < lamda2) .and. (lamda2 < eps*ctil) )
            lamda2 = HALF*(eps*ctil + lamda2*lamda2/(eps*ctil))
        else where
            lamda2 = sqrt(lamda2*lamda2)
        end where

        where ( (-eps*ctil < lamda3) .and. (lamda3 < eps*ctil) )
            lamda3 = HALF*(eps*ctil + lamda3*lamda3/(eps*ctil))
        else where
            lamda3 = sqrt(lamda3*lamda3)
        end where





        C1   = abs(lamda1)*( delp - rtil*ctil*delvmag)/(TWO*(ctil2))
        C2_a = abs(lamda2)*(delr - delp/(ctil2))
        C2_b = abs(lamda2)*rtil
        C3   = abs(lamda3)*( delp + rtil*ctil*delvmag)/(TWO*(ctil2))

        integrand = delr
        integrand = ZERO


        !================================
        !       MASS FLUX
        !================================
        upwind = C1 + C2_a + C3

        integrand = HALF*(upwind*normx*unormx + upwind*normy*unormy + upwind*normz*unormz)


        call worker%integrate_boundary('Density',integrand)


        !================================
        !       X-MOMENTUM FLUX
        !================================
        upwind = C1*(util - ctil*unormx)  +  C2_a*util  +  C2_b*(delu - delvmag*unormx)  +  C3*(util + ctil*unormx)

        integrand = HALF*(upwind*normx*unormx + upwind*normy*unormy + upwind*normz*unormz)


        call worker%integrate_boundary('Momentum-1',integrand)


        !================================
        !       Y-MOMENTUM FLUX
        !================================
        upwind = C1*(vtil - ctil*unormy)  +  C2_a*vtil  +  C2_b*(delv - delvmag*unormy)  +  C3*(vtil + ctil*unormy)

        integrand = HALF*(upwind*normx*unormx + upwind*normy*unormy + upwind*normz*unormz)


        call worker%integrate_boundary('Momentum-2',integrand)

        !================================
        !       Z-MOMENTUM FLUX
        !================================
        upwind = C1*(wtil - ctil*unormz)  +  C2_a*wtil  +  C2_b*(delw - delvmag*unormz)  +  C3*(wtil + ctil*unormz)

        integrand = HALF*(upwind*normx*unormx + upwind*normy*unormy + upwind*normz*unormz)


        call worker%integrate_boundary('Momentum-3',integrand)

        !================================
        !          ENERGY FLUX
        !================================
        upwind = C1*(Htil - ctil*vmagtil)  +  C2_a*(qtil2/TWO)  +  C2_b*(util*delu + vtil*delv + wtil*delw - vmagtil*delvmag)  +  C3*(Htil + ctil*vmagtil)

        integrand = HALF*(upwind*normx*unormx + upwind*normy*unormy + upwind*normz*unormz)


        call worker%integrate_boundary('Energy',integrand)


    end subroutine compute
    !**********************************************************************************************













end module euler_ale_roe_operator
