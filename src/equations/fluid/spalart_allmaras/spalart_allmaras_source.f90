module spalart_allmaras_source
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,SIX,HALF,PI
    use mod_spalart_allmaras,   only: SA_c_b1, SA_c_b2, SA_kappa, SA_sigma,         &
                                      SA_c_w1, SA_c_w2, SA_c_w3, SA_c_v1, SA_c_v2,  &
                                      SA_c_v3, SA_c_t3, SA_c_t4, SA_c_n1,           &
                                      SA_Pr_t, SA_rlim, SA_b    

    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    implicit none

    private

    
    !> Volume flux for Fluid Viscous Terms.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/23/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: spalart_allmaras_source_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type spalart_allmaras_source_operator_t
    !******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/23/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(spalart_allmaras_source_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("Fluid Spalart-Allmaras Source Operator")

        ! Set operator type
        call self%set_operator_type("Volume Diffusive Operator")

        ! Set operator equations
        call self%set_equation("Density * NuTilde")

    end subroutine init
    !********************************************************************************



    !> Volume flux routine for Fluid Viscous Terms.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/23/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(spalart_allmaras_source_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop

        ! Equation indices
        integer(ik)    :: irho, irhou, irhov, irhow, irhoE, irho_nutilde


        type(AD_D), allocatable, dimension(:) ::                                &
            rho, rhou, rhov, rhow, rhoE, rho_nutilde, p, T, u, v, w, invrho,    &
            gam, mu, nu, lamda,                                                 &
            drho_dx, drhou_dx, drhov_dx, drhow_dx, drhoE_dx,                    &
            drho_dy, drhou_dy, drhov_dy, drhow_dy, drhoE_dy,                    &
            drho_dz, drhou_dz, drhov_dz, drhow_dz, drhoE_dz,                    &
            drho_nutilde_dx, drho_nutilde_dy, drho_nutilde_dz,                  &
            du_dx,   dv_dx,    dw_dx,    dT_dx,                                 &
            du_dy,   dv_dy,    dw_dy,    dT_dy,                                 &
            du_dz,   dv_dz,    dw_dz,    dT_dz,                                 &
            du_drho, du_drhou, dv_drho,  dv_drhov, dw_drho, dw_drhow,           &
            vorticity2, vorticity, vorticity_bar, vorticity_mod,                &
            f_v1, f_v2, f_t2, f_w, f_n1, production, destruction,               &
            dnutilde_drho, dnutilde_drho_nutilde, chi, r, rbar, g, mu_t,        &
            nutilde, dnutilde_dx, dnutilde_dy, dnutilde_dz,                     &
            source

        real(rk), allocatable, dimension(:) :: dwall

        real(rk)    :: const, epsilon_vorticity

        !
        ! Get equation indices
        !
        irho         = prop%get_equation_index("Density"   )
        irhou        = prop%get_equation_index("X-Momentum")
        irhov        = prop%get_equation_index("Y-Momentum")
        irhow        = prop%get_equation_index("Z-Momentum")
        irhoE        = prop%get_equation_index("Energy"    )
        irho_nutilde = prop%get_equation_index("Density * NuTilde")



        !
        ! Interpolate solution to quadrature nodes
        !
        rho         = worker%get_element_variable(irho,        'value')
        rhou        = worker%get_element_variable(irhou,       'value')
        rhov        = worker%get_element_variable(irhov,       'value')
        rhow        = worker%get_element_variable(irhow,       'value')
        rhoE        = worker%get_element_variable(irhoE,       'value')
        rho_nutilde = worker%get_element_variable(irho_nutilde,'value')






        !
        ! Interpolate solution gradients to quadrature nodes
        !
        drho_dx  = worker%get_element_variable(irho,  'ddx+lift')
        drho_dy  = worker%get_element_variable(irho,  'ddy+lift')
        drho_dz  = worker%get_element_variable(irho,  'ddz+lift')

        drhou_dx = worker%get_element_variable(irhou, 'ddx+lift')
        drhou_dy = worker%get_element_variable(irhou, 'ddy+lift')
        drhou_dz = worker%get_element_variable(irhou, 'ddz+lift')

        drhov_dx = worker%get_element_variable(irhov, 'ddx+lift')
        drhov_dy = worker%get_element_variable(irhov, 'ddy+lift')
        drhov_dz = worker%get_element_variable(irhov, 'ddz+lift')

        drhow_dx = worker%get_element_variable(irhow, 'ddx+lift')
        drhow_dy = worker%get_element_variable(irhow, 'ddy+lift')
        drhow_dz = worker%get_element_variable(irhow, 'ddz+lift')

        drhoE_dx = worker%get_element_variable(irhoE, 'ddx+lift')
        drhoE_dy = worker%get_element_variable(irhoE, 'ddy+lift')
        drhoE_dz = worker%get_element_variable(irhoE, 'ddz+lift')

        drho_nutilde_dx = worker%get_element_variable(irho_nutilde, 'ddx+lift')
        drho_nutilde_dy = worker%get_element_variable(irho_nutilde, 'ddy+lift')
        drho_nutilde_dz = worker%get_element_variable(irho_nutilde, 'ddz+lift')



        !
        ! Divide by density
        !
        invrho  = ONE/rho
        u       = rhou*invrho
        v       = rhov*invrho
        w       = rhow*invrho
        nutilde = rho_nutilde*invrho



        !
        ! Compute model values
        !
        p   = prop%fluid%compute_pressure(rho,rhou,rhov,rhow,rhoE)
        gam = prop%fluid%compute_gamma(rho,rhou,rhov,rhow,rhoE)
        T   = prop%fluid%compute_temperature(rho,rhou,rhov,rhow,rhoE)
        mu  = prop%fluid%compute_viscosity_dynamic(T)
        nu  = mu*invrho



        !
        ! Compute turbulence viscosity
        !
        chi  = nutilde/nu
        f_v1 = (chi*chi*chi)/(chi*chi*chi + SA_c_v1*SA_c_v1*SA_c_v1)


        mu_t = rho_nutilde
        mu_t = ZERO
        where(nutilde >= ZERO)
            mu_t = rho_nutilde * f_v1
        end where
        


        !
        ! Compute f_n1
        !
        f_n1 = rho_nutilde
        f_n1 = ONE
        where(nutilde < ZERO)
            f_n1 = (SA_c_n1 + chi*chi*chi)/(SA_c_n1 - chi*chi*chi)
        end where



        !
        ! Compute velocity jacobians
        !
        du_drho  = -invrho*invrho*rhou
        du_drhou =  invrho

        dv_drho  = -invrho*invrho*rhov
        dv_drhov =  invrho

        dw_drho  = -invrho*invrho*rhow
        dw_drhow =  invrho




        !
        ! Compute velocity gradients
        !
        du_dx = du_drho*drho_dx  +  du_drhou*drhou_dx
        du_dy = du_drho*drho_dy  +  du_drhou*drhou_dy
        du_dz = du_drho*drho_dz  +  du_drhou*drhou_dz

        dv_dx = dv_drho*drho_dx  +  dv_drhov*drhov_dx
        dv_dy = dv_drho*drho_dy  +  dv_drhov*drhov_dy
        dv_dz = dv_drho*drho_dz  +  dv_drhov*drhov_dz

        dw_dx = dw_drho*drho_dx  +  dw_drhow*drhow_dx
        dw_dy = dw_drho*drho_dy  +  dw_drhow*drhow_dy
        dw_dz = dw_drho*drho_dz  +  dw_drhow*drhow_dz



        !
        ! Compute vorticity and modified vorticity
        !
        vorticity2 =  (dw_dy - dv_dz)**TWO  +  (du_dz - dw_dx)**TWO  +  (dv_dx - du_dy)**TWO 
        
        epsilon_vorticity = 1.e-6_rk
        where(vorticity2 < epsilon_vorticity)
            vorticity = HALF*(epsilon_vorticity + vorticity2/epsilon_vorticity)
        else where
            vorticity = sqrt(vorticity2)
        end where


        f_v2 = ONE - (chi/(ONE+chi*f_v1))
        vorticity_bar = (nutilde/(SA_kappa*SA_kappa*dwall*dwall))*f_v2


        where (vorticity_bar >= -SA_c_v2*vorticity)
            vorticity_mod = vorticity + vorticity_bar
        else where
            vorticity_mod = vorticity + vorticity*(SA_c_v2*SA_c_v2*vorticity + SA_c_v3*vorticity_bar)/( (SA_c_v3 - TWO*SA_c_v2)*vorticity - vorticity_bar ) 
        end where




        !
        ! Compute f_t2, f_w, g, r
        !
        rbar = nutilde/(vorticity_mod * SA_kappa * SA_kappa * dwall * dwall )

        !r = min(SA_rlim,rbar)
        r = SA_rlim - (SA_rlim - rbar)*(atan(SA_b * (SA_rlim - rbar))/PI  + HALF)  + atan(SA_b)/PI  -  HALF

        g = r + SA_c_w2*(r**SIX  -  r)

        f_w = g*(((ONE + SA_c_w3**SIX)/(g**SIX + SA_c_w3**SIX))**(ONE/SIX))

        f_t2 = SA_c_t3*exp(-SA_c_t4*chi*chi)



        !
        ! Compute Production, Destruction
        !
        where ( nutilde >= ZERO )
            production = SA_c_b1*(ONE - f_t2)*vorticity_mod*nutilde
        else where
            production = SA_C_b1*(ONE - SA_c_t3)*vorticity*nutilde
        end where


        where ( nutilde >= ZERO )
            destruction = (SA_c_w1*f_w - (SA_c_b1/(SA_kappa*SA_kappa))*f_t2) * (nutilde/dwall)**TWO
        else where
            destruction = -SA_c_w1 * (nutilde/dwall)**TWO
        end where


        !
        ! Compute partial derivatives of nutilde
        !
        dnutilde_drho         = -invrho*invrho*rho_nutilde
        dnutilde_drho_nutilde =  invrho


        !
        ! Compute cartesian derivatives
        !
        dnutilde_dx = dnutilde_drho*drho_dx  +  dnutilde_drho_nutilde*drho_nutilde_dx
        dnutilde_dy = dnutilde_drho*drho_dy  +  dnutilde_drho_nutilde*drho_nutilde_dy
        dnutilde_dz = dnutilde_drho*drho_dz  +  dnutilde_drho_nutilde*drho_nutilde_dz



        !========================================================================
        !                       Spalart-Allmaras Source Term
        !========================================================================
        source = -rho*(production-destruction)  -  &
                  (SA_c_b2/SA_sigma)*rho*(dnutilde_dx*dnutilde_dx + dnutilde_dy*dnutilde_dy + dnutilde_dz*dnutilde_dz)  +  &
                  (ONE/SA_sigma)*(nu + f_n1*nutilde)*(drho_dx*dnutilde_dx + drho_dy*dnutilde_dy + drho_dz*dnutilde_dz)

        call worker%integrate_volume(irho_nutilde, source)

    end subroutine compute
    !*********************************************************************************************************






end module spalart_allmaras_source
