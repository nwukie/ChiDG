module euler_roe_operator
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF
    use mod_fluid,              only: gam
    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
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
    type, extends(operator_t), public :: euler_roe_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type euler_roe_operator_t
    !*******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_roe_operator_t),   intent(inout)    :: self

        ! Set operator name
        call self%set_name("Euler Roe Flux")

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
        class(euler_roe_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                  &
            density_m,      density_p,                              &
            mom1_m,         mom1_p,                                 &
            mom2_m,         mom2_p,                                 &
            mom3_m,         mom3_p,                                 &
            energy_m,       energy_p,                               &
            enthalpy_m,     enthalpy_p,                             &
            invdensity_m,   invdensity_p,                           &
            p_m,            p_p,                                    &
            un_m,           un_p,                                   &
            rtil, util, vtil, wtil, vmagtil, Htil, ctil, qtil2,     &
            C1,  C2_a, C2_b,  C3,                                   &
            u_m, v_m, w_m,                                          &
            u_p, v_p, w_p,                                          &
            vmag_p, vmag_m,                                         &
            delr,   delp,   delvmag, delu, delv, delw,              &
            lamda1, lamda2, lamda3,                                 &
            sqrt_rhom, sqrt_rhop, sqrt_rhom_plus_rhop, ctil2, upwind

        real(rk), allocatable, dimension(:) ::              &
            norm_1,         norm_2,         norm_3,         &
            unorm_1,        unorm_2,        unorm_3,        &
            unorm_1_ale,    unorm_2_ale,    unorm_3_ale,    &
            ale_area_ratio, r, grid_vel_n

        real(rk), allocatable, dimension(:,:) :: grid_vel

        real(rk) :: eps


        grid_vel = worker%get_grid_velocity_face('face interior')
        ale_area_ratio = worker%get_area_ratio()


        !
        ! Interpolate solution to quadrature nodes
        !
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


        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','face interior') 
            mom2_m = mom2_m / r
            mom2_p = mom2_p / r
        end if



        norm_1  = worker%normal(1)
        norm_2  = worker%normal(2)
        norm_3  = worker%normal(3)

        unorm_1 = worker%unit_normal(1)
        unorm_2 = worker%unit_normal(2)
        unorm_3 = worker%unit_normal(3)

        unorm_1_ale = worker%unit_normal_ale(1)
        unorm_2_ale = worker%unit_normal_ale(2)
        unorm_3_ale = worker%unit_normal_ale(3)




        !
        ! Compute pressure and gamma
        !
        p_m = worker%get_field('Pressure', 'value', 'face interior')
        p_p = worker%get_field('Pressure', 'value', 'face exterior')

        invdensity_m = ONE/density_m
        invdensity_p = ONE/density_p

        !
        ! Compute enthalpy
        !
        enthalpy_m = (energy_m + p_m)*invdensity_m
        enthalpy_p = (energy_p + p_p)*invdensity_p


        !
        ! Compute velocity components
        !
        u_m = mom1_m*invdensity_m 
        v_m = mom2_m*invdensity_m 
        w_m = mom3_m*invdensity_m 
        vmag_m = u_m*unorm_1_ale + v_m*unorm_2_ale + w_m*unorm_3_ale

        u_p = mom1_p*invdensity_p 
        v_p = mom2_p*invdensity_p 
        w_p = mom3_p*invdensity_p 
        vmag_p = u_p*unorm_1_ale + v_p*unorm_2_ale + w_p*unorm_3_ale


        !
        ! Compute Roe-averaged variables
        !
        sqrt_rhom = sqrt(density_m)
        sqrt_rhop = sqrt(density_p)
        sqrt_rhom_plus_rhop = sqrt_rhom + sqrt_rhop
        rtil =  sqrt(density_p * density_m)                                             ! Roe-averaged density
        util = (sqrt_rhom*u_m        + sqrt_rhop*u_p       ) / (sqrt_rhom_plus_rhop)    ! Roe-averaged u-velocity
        vtil = (sqrt_rhom*v_m        + sqrt_rhop*v_p       ) / (sqrt_rhom_plus_rhop)    ! Roe-averaged v-velocity
        wtil = (sqrt_rhom*w_m        + sqrt_rhop*w_p       ) / (sqrt_rhom_plus_rhop)    ! Roe-averaged w-velocity
        Htil = (sqrt_rhom*enthalpy_m + sqrt_rhop*enthalpy_p) / (sqrt_rhom_plus_rhop)    ! Roe-averaged Enthalpy

        vmagtil = util*unorm_1_ale + vtil*unorm_2_ale + wtil*unorm_3_ale   ! Magnitude of Roe-averaged velocity in the face normal direction
        qtil2   = util**TWO + vtil**TWO + wtil**TWO
        ctil = sqrt((gam - ONE)*(Htil - HALF*qtil2))                   ! Roe-averaged speed of sound
        ctil2 = ctil**TWO



        !
        ! Compute jump terms
        !
        delr    = (density_m - density_p)
        delu    = (u_m - u_p)
        delv    = (v_m - v_p)
        delw    = (w_m - w_p)
        delvmag = (vmag_m - vmag_p)
        delp    = (p_m - p_p)


        !
        ! Limit wave speeds for entropy fix
        !
        !lamda1 = abs(vmagtil - ctil) + sqrt(grid_vel(:,1)**TWO+grid_vel(:,2)**TWO+grid_vel(:,3)**TWO)
        !lamda2 = abs(vmagtil)        + sqrt(grid_vel(:,1)**TWO+grid_vel(:,2)**TWO+grid_vel(:,3)**TWO)
        !lamda3 = abs(vmagtil + ctil) + sqrt(grid_vel(:,1)**TWO+grid_vel(:,2)**TWO+grid_vel(:,3)**TWO)

        grid_vel_n = grid_vel(:,1)*unorm_1_ale  +  grid_vel(:,2)*unorm_2_ale  +  grid_vel(:,3)*unorm_3_ale
        lamda1 = abs(vmagtil + grid_vel_n - ctil)
        lamda2 = abs(vmagtil + grid_vel_n)
        lamda3 = abs(vmagtil + grid_vel_n + ctil)

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



        !================================
        !       Mass flux
        !================================
        upwind = HALF*(C1 + C2_a + C3)

        call worker%integrate_boundary_upwind('Density',upwind)


        !================================
        !       Momentum-1 flux
        !================================
        upwind = HALF*(C1*(util - ctil*unorm_1_ale)  +  C2_a*util  +  C2_b*(delu - delvmag*unorm_1_ale)  +  C3*(util + ctil*unorm_1_ale))

        call worker%integrate_boundary_upwind('Momentum-1',upwind)


        !================================
        !       Momentum-2 flux
        !================================
        upwind = HALF*(C1*(vtil - ctil*unorm_2_ale)  +  C2_a*vtil  +  C2_b*(delv - delvmag*unorm_2_ale)  +  C3*(vtil + ctil*unorm_2_ale))

        ! Convert to tangential to angular momentum flux
        if (worker%coordinate_system() == 'Cylindrical') then
            upwind = upwind * r
        end if

        call worker%integrate_boundary_upwind('Momentum-2',upwind)

        !================================
        !       Momentum-3 flux
        !================================
        upwind = HALF*(C1*(wtil - ctil*unorm_3_ale)  +  C2_a*wtil  +  C2_b*(delw - delvmag*unorm_3_ale)  +  C3*(wtil + ctil*unorm_3_ale))

        call worker%integrate_boundary_upwind('Momentum-3',upwind)

        !================================
        !       Energy flux
        !================================
        upwind = HALF*(C1*(Htil - ctil*vmagtil)  +  C2_a*(qtil2/TWO)  +  C2_b*(util*delu + vtil*delv + wtil*delw - vmagtil*delvmag)  +  C3*(Htil + ctil*vmagtil))

        call worker%integrate_boundary_upwind('Energy',upwind)


    end subroutine compute
    !**********************************************************************************************













end module euler_roe_operator
