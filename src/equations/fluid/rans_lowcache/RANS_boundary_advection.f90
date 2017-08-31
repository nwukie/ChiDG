module RANS_boundary_advection
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF
    use mod_fluid,              only: omega, gam
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
    type, extends(operator_t), public :: RANS_boundary_advection_t

!        real(rk)    :: gam = 1.4_rk
!        real(rk)    :: R   = 287.15_rk
!        real(rk)    :: Cp  = 1003.0_rk
!        real(rk)    :: Pr  = 0.72_rk

    contains

        procedure   :: init
        procedure   :: compute

    end type RANS_boundary_advection_t
    !*******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(RANS_boundary_advection_t),   intent(inout)    :: self

        ! Set operator name
        call self%set_name("RANS Boundary Advection")

        ! Set operator type
        call self%set_operator_type("Boundary Advective Flux")

        ! Set operator equations
        call self%add_primary_field('Density'   )
        call self%add_primary_field('Momentum-1')
        call self%add_primary_field('Momentum-2')
        call self%add_primary_field('Momentum-3')
        call self%add_primary_field('Energy'    )
        call self%add_primary_field('Density * NuTilde')

    end subroutine init
    !********************************************************************************




    !> Compute Roe approximate Riemann upwind flux
    !! 
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!---------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(RANS_boundary_advection_t),    intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                  &
            density_m,      density_p,                              &
            mom1_m,         mom1_p,                                 &
            mom2_m,         mom2_p,                                 &
            mom3_m,         mom3_p,                                 &
            energy_m,       energy_p,                               &
            density_nutilde_m, density_nutilde_p,                   &
            enthalpy_m,     enthalpy_p,                             &
            invdensity_m,   invdensity_p,                           &
            diss_m,         diss_p,                                 &
            p_m,            p_p,                                    &
            un_m,           un_p,                                   &
            a_m,            a_p,                                    &
            u_m, v_m, w_m,                                    &
            u_p, v_p, w_p,                                    &
            rtil, util, vtil, wtil, vmagtil, Htil, ctil, qtil2,     &
            vtil_t, v_m_t, v_p_t, vmagtil_t,                        &
            integrand,  upwind,     wave,                           &
            C1,  C2_a, C2_b,  C3,                                   &
            u_m, v_m, w_m,                                          &
            u_p, v_p, w_p,                                          &
            c_m, c_p, t_m, t_p, diff,                               &
            vmag_p, vmag_m,                                         &
            delr,   delp,   delvmag, delu, delv, delw,              &
            lamda1, lamda2, lamda3,                                 &
            sqrt_rhom, sqrt_rhop, sqrt_rhom_plus_rhop, ctil2,       &
            flux_avg_1, flux_avg_2, flux_avg_3,                     &
            flux_1_m, flux_2_m, flux_3_m,                           &
            flux_1_p, flux_2_p, flux_3_p,                           &
            flux_1,   flux_2,   flux_3

        real(rk), allocatable, dimension(:) :: &
            norm_1, norm_2, norm_3, unorm_1, unorm_2, unorm_3, r, area

        real(rk) :: eps



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

        density_nutilde_m = worker%get_primary_field_face('Density * NuTilde', 'value', 'face interior')
        density_nutilde_p = worker%get_primary_field_face('Density * NuTilde', 'value', 'face exterior')


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
        ! Get scaled and unit normal vectors
        !
        norm_1  = worker%normal(1)
        norm_2  = worker%normal(2)
        norm_3  = worker%normal(3)

        unorm_1 = worker%unit_normal(1)
        unorm_2 = worker%unit_normal(2)
        unorm_3 = worker%unit_normal(3)



        !
        ! Compute pressure and gamma
        !
        invdensity_m = ONE/density_m
        invdensity_p = ONE/density_p
        p_m = (self%gam-ONE)*(energy_m - HALF*( (mom1_m*mom1_m) + (mom2_m*mom2_m) + (mom3_m*mom3_m) )*invdensity_m )
        p_p = (self%gam-ONE)*(energy_p - HALF*( (mom1_p*mom1_p) + (mom2_p*mom2_p) + (mom3_p*mom3_p) )*invdensity_p )
        t_m = p_m/(density_m*self%R)
        t_p = p_p/(density_p*self%R)
        enthalpy_m = (energy_m + p_m)*invdensity_m
        enthalpy_p = (energy_p + p_p)*invdensity_p


        !
        ! Compute velocity components
        !
        u_m = mom1_m*invdensity_m
        v_m = mom2_m*invdensity_m
        w_m = mom3_m*invdensity_m

        u_p = mom1_p*invdensity_p
        v_p = mom2_p*invdensity_p
        w_p = mom3_p*invdensity_p


        !
        ! Compute transport velocities 
        !
        v_m_t = v_m - omega*r
        v_p_t = v_p - omega*r


        vmag_m = u_m*unorm_1 + v_m*unorm_2 + w_m*unorm_3
        vmag_p = u_p*unorm_1 + v_p*unorm_2 + w_p*unorm_3


        !
        ! Compute Roe-averaged variables
        !
        sqrt_rhom = sqrt(density_m)
        sqrt_rhop = sqrt(density_p)
        sqrt_rhom_plus_rhop = sqrt_rhom + sqrt_rhop
        rtil =  sqrt(density_p * density_m)                                     ! Roe-averaged density
        util = (sqrt_rhom*u_m + sqrt_rhop*u_p) / (sqrt_rhom_plus_rhop)          ! Roe-averaged u-velocity
        vtil = (sqrt_rhom*v_m + sqrt_rhop*v_p) / (sqrt_rhom_plus_rhop)          ! Roe-averaged v-velocity
        vtil_t = (sqrt_rhom*v_m_t + sqrt_rhop*v_p_t) / (sqrt_rhom_plus_rhop)    ! Roe-averaged v-velocity
        wtil = (sqrt_rhom*w_m + sqrt_rhop*w_p) / (sqrt_rhom_plus_rhop)          ! Roe-averaged w-velocity
        Htil = (sqrt_rhom*enthalpy_m + sqrt_rhop*enthalpy_p) / (sqrt_rhom_plus_rhop)    ! Roe-averaged Enthalpy

        ! Magnitude of Roe-averaged velocity in the face-normal direction
        vmagtil   = util*unorm_1 + vtil  *unorm_2 + wtil*unorm_3
        vmagtil_t = util*unorm_1 + vtil_t*unorm_2 + wtil*unorm_3
        qtil2     = util**TWO + vtil**TWO + wtil**TWO

        !& HARDCODED GAMMA
        ctil = sqrt((1.4_rk - ONE)*(Htil - HALF*qtil2))                   ! Roe-averaged speed of sound
        ctil2 = ctil**TWO



        !
        ! Compute jump terms
        !
        delr    = (density_m - density_p)
        delu    = (u_m       - u_p      )
        delv    = (v_m       - v_p      )
        delw    = (w_m       - w_p      )
        delvmag = (vmag_m    - vmag_p   )
        delp    = (p_m       - p_p      )


        !
        ! Limit wave speeds for entropy fix
        !
        lamda1 = abs(vmagtil_t - ctil)
        lamda2 = abs(vmagtil_t       )
        lamda3 = abs(vmagtil_t + ctil)

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



        !
        ! Get fluid advection velocity
        !
        u_m = worker%get_model_field_face('Advection Velocity-1', 'value', 'face interior')
        v_m = worker%get_model_field_face('Advection Velocity-2', 'value', 'face interior')
        w_m = worker%get_model_field_face('Advection Velocity-3', 'value', 'face interior')

        u_p = worker%get_model_field_face('Advection Velocity-1', 'value', 'face exterior')
        v_p = worker%get_model_field_face('Advection Velocity-2', 'value', 'face exterior')
        w_p = worker%get_model_field_face('Advection Velocity-3', 'value', 'face exterior')









        !
        ! Compute differential areas
        !
        area = sqrt(norm_1**TWO + norm_2**TWO + norm_3**TWO)


        !=================================================
        ! mass flux
        !=================================================
        flux_1_m = (density_m * u_m)
        flux_2_m = (density_m * v_m)
        flux_3_m = (density_m * w_m)

        flux_1_p = (density_p * u_p)
        flux_2_p = (density_p * v_p)
        flux_3_p = (density_p * w_p)

        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)


        upwind = HALF*area*(C1 + C2_a + C3)

        !integrand = HALF*(upwind*norm_1*unorm_1 + upwind*norm_2*unorm_2 + upwind*norm_3*unorm_3)
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3) + upwind


        call worker%integrate_boundary('Density',integrand)


        !=================================================
        ! momentum-1 flux
        !=================================================
        flux_1_m = (mom1_m * u_m) + p_m
        flux_2_m = (mom1_m * v_m)
        flux_3_m = (mom1_m * w_m)

        flux_1_p = (mom1_p * u_p) + p_p
        flux_2_p = (mom1_p * v_p)
        flux_3_p = (mom1_p * w_p)

        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)

        upwind = HALF*area*(C1*(util - ctil*unorm_1)  +  C2_a*util  +  C2_b*(delu - delvmag*unorm_1)  +  C3*(util + ctil*unorm_1))

        !integrand = HALF*(upwind*norm_1*unorm_1 + upwind*norm_2*unorm_2 + upwind*norm_3*unorm_3)
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)  +  upwind

        call worker%integrate_boundary('Momentum-1',integrand)


        !=================================================
        ! momentum-2 flux
        !=================================================
        flux_1_m = (mom2_m * u_m)
        flux_2_m = (mom2_m * v_m) + p_m
        flux_3_m = (mom2_m * w_m)

        flux_1_p = (mom2_p * u_p)
        flux_2_p = (mom2_p * v_p) + p_p
        flux_3_p = (mom2_p * w_p)

        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)

        upwind = HALF*area*(C1*(vtil - ctil*unorm_2)  +  C2_a*vtil  +  C2_b*(delv - delvmag*unorm_2)  +  C3*(vtil + ctil*unorm_2))

        !integrand = HALF*(upwind*norm_1*unorm_1 + upwind*norm_2*unorm_2 + upwind*norm_3*unorm_3)
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)  +  upwind

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
        flux_1_m = (mom3_m * u_m)
        flux_2_m = (mom3_m * v_m)
        flux_3_m = (mom3_m * w_m) + p_m

        flux_1_p = (mom3_p * u_p)
        flux_2_p = (mom3_p * v_p)
        flux_3_p = (mom3_p * w_p) + p_p


        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)


        upwind = HALF*area*(C1*(wtil - ctil*unorm_3)  +  C2_a*wtil  +  C2_b*(delw - delvmag*unorm_3)  +  C3*(wtil + ctil*unorm_3))

        !integrand = HALF*(upwind*norm_1*unorm_1 + upwind*norm_2*unorm_2 + upwind*norm_3*unorm_3)
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)  +  upwind

        call worker%integrate_boundary('Momentum-3',integrand)

        !=================================================
        ! energy flux
        !=================================================
        flux_1_m = (density_m * enthalpy_m * u_m)
        flux_2_m = (density_m * enthalpy_m * v_m)  +  omega*r*p_m
        flux_3_m = (density_m * enthalpy_m * w_m)

        flux_1_p = (density_p * enthalpy_p * u_p)
        flux_2_p = (density_p * enthalpy_p * v_p)  +  omega*r*p_p
        flux_3_p = (density_p * enthalpy_p * w_p)
        
        flux_1 = (flux_1_m + flux_1_p)
        flux_2 = (flux_2_m + flux_2_p)
        flux_3 = (flux_3_m + flux_3_p)



        upwind = HALF*area*(C1*(Htil - ctil*vmagtil)  +  C2_a*(qtil2/TWO)  +  C2_b*(util*delu + vtil*delv + wtil*delw - vmagtil*delvmag)  +  C3*(Htil + ctil*vmagtil))

        !integrand = HALF*(upwind*norm_1*unorm_1 + upwind*norm_2*unorm_2 + upwind*norm_3*unorm_3)
        integrand = HALF*(flux_1*norm_1 + flux_2*norm_2 + flux_3*norm_3)  +  upwind

        call worker%integrate_boundary('Energy',integrand)



        !=================================================
        ! turbulence flux
        !=================================================


        !
        ! Compute average flux
        ! 
        flux_avg_1 = HALF*(density_nutilde_m*u_m  +  density_nutilde_p*u_p)
        flux_avg_2 = HALF*(density_nutilde_m*v_m  +  density_nutilde_p*v_p)
        flux_avg_3 = HALF*(density_nutilde_m*w_m  +  density_nutilde_p*w_p)

        !
        ! Compute maximum wave speed
        !
        un_m = u_m*unorm_1 + v_m*unorm_2 + w_m*unorm_3
        un_p = u_p*unorm_1 + v_p*unorm_2 + w_p*unorm_3


        c_m = sqrt(self%gam * self%R * t_m)
        c_p = sqrt(self%gam * self%R * t_p)

        diss_m = abs(un_m) + c_m
        diss_p = abs(un_p) + c_p


        !
        ! Compute Lax-Friedrichs upwind flux
        !
        diff   = (density_nutilde_m - density_nutilde_p)
        !flux_1 = flux_avg_1 + max(abs(diss_m),abs(diss_p))*HALF*diff
        !flux_2 = flux_avg_2 + max(abs(diss_m),abs(diss_p))*HALF*diff
        !flux_3 = flux_avg_3 + max(abs(diss_m),abs(diss_p))*HALF*diff


        integrand = flux_avg_1*norm_1 + flux_avg_2*norm_2 + flux_avg_3*norm_3

        integrand = integrand + max(abs(diss_m),abs(diss_p))*HALF*diff*norm_1*unorm_1
        integrand = integrand + max(abs(diss_m),abs(diss_p))*HALF*diff*norm_2*unorm_2
        integrand = integrand + max(abs(diss_m),abs(diss_p))*HALF*diff*norm_3*unorm_3

        !
        ! Integrate flux
        !
        call worker%integrate_boundary('Density * NuTilde',integrand)




    end subroutine compute
    !**********************************************************************************************













end module RANS_boundary_advection
