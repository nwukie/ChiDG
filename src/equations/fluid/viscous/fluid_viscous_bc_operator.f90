module fluid_viscous_bc_operator
    use mod_kinds,          only: ik, rk
    use mod_constants,      only: BC, HALF, ONE, TWO
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !------------------------------------------------------------------------------------
    type, public, extends(operator_t)   :: fluid_viscous_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type fluid_viscous_bc_operator_t
    !*************************************************************************************






contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(fluid_viscous_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Fluid Viscous BC Operator")

        !
        ! Set operator type
        !
        call self%set_operator_type("BC Diffusive Flux")

        !
        ! Set operator equations
        !
        call self%set_equation("Density"   )
        call self%set_equation("X-Momentum")
        call self%set_equation("Y-Momentum")
        call self%set_equation("Z-Momentum")
        call self%set_equation("Energy"    )

    end subroutine init
    !********************************************************************************





    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!  @param[in]      face    face_info_t containing indices on location and misc information
    !!  @param[in]      flux    function_into_t containing info on the function being computed
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(fluid_viscous_bc_operator_t), intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::                                &
            rho, rhou, rhov, rhow, rhoE, p, T, u, v, w, invrho, gam, mu, lamda, &
            drho_dx, drhou_dx, drhov_dx, drhow_dx, drhoE_dx,                    &
            drho_dy, drhou_dy, drhov_dy, drhow_dy, drhoE_dy,                    &
            drho_dz, drhou_dz, drhov_dz, drhow_dz, drhoE_dz,                    &
            du_dx,   dv_dx,    dw_dx,    dT_dx,                                 &
            du_dy,   dv_dy,    dw_dy,    dT_dy,                                 &
            du_dz,   dv_dz,    dw_dz,    dT_dz,                                 &
            du_drho, du_drhou, dv_drho,  dv_drhov, dw_drho, dw_drhow,           &
            dT_drho, dT_drhou, dT_drhov, dT_drhow, dT_drhoE,                    &
            dp_drho, dp_drhou, dp_drhov, dp_drhow, dp_drhoE,                    &
            dke_drho, dke_drhou, dke_drhov, dke_drhow,                          &
            tau_xx, tau_yy, tau_zz, tau_xy, tau_xz, tau_yz,                     &
            flux_x, flux_y, flux_z, integrand

        real(rk),   allocatable, dimension(:)   ::          &
            normx, normy, normz

        real(rk) :: const

        !
        ! Get equation indices
        !
        irho  = prop%get_equation_index("Density"   )
        irhou = prop%get_equation_index("X-Momentum")
        irhov = prop%get_equation_index("Y-Momentum")
        irhow = prop%get_equation_index("Z-Momentum")
        irhoE = prop%get_equation_index("Energy"    )



        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        rho  = worker%get_face_variable(irho,  'value', BC)
        rhou = worker%get_face_variable(irhou, 'value', BC)
        rhov = worker%get_face_variable(irhov, 'value', BC)
        rhow = worker%get_face_variable(irhow, 'value', BC)
        rhoE = worker%get_face_variable(irhoE, 'value', BC)


        !
        ! Compute model values
        !
        p   = prop%fluid%compute_pressure(rho,rhou,rhov,rhow,rhoE)
        gam = prop%fluid%compute_gamma(rho,rhou,rhov,rhow,rhoE)


        !
        ! Interpolate solution gradients to quadrature nodes
        !
        drho_dx  = worker%get_face_variable(irho,  'ddx+lift', BC)
        drho_dy  = worker%get_face_variable(irho,  'ddy+lift', BC)
        drho_dz  = worker%get_face_variable(irho,  'ddz+lift', BC)

        drhou_dx = worker%get_face_variable(irhou, 'ddx+lift', BC)
        drhou_dy = worker%get_face_variable(irhou, 'ddy+lift', BC)
        drhou_dz = worker%get_face_variable(irhou, 'ddz+lift', BC)

        drhov_dx = worker%get_face_variable(irhov, 'ddx+lift', BC)
        drhov_dy = worker%get_face_variable(irhov, 'ddy+lift', BC)
        drhov_dz = worker%get_face_variable(irhov, 'ddz+lift', BC)

        drhow_dx = worker%get_face_variable(irhow, 'ddx+lift', BC)
        drhow_dy = worker%get_face_variable(irhow, 'ddy+lift', BC)
        drhow_dz = worker%get_face_variable(irhow, 'ddz+lift', BC)

        drhoE_dx = worker%get_face_variable(irhoE, 'ddx+lift', BC)
        drhoE_dy = worker%get_face_variable(irhoE, 'ddy+lift', BC)
        drhoE_dz = worker%get_face_variable(irhoE, 'ddz+lift', BC)


        invrho = ONE/rho




        !
        ! Get normal vector
        !
        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)




        !
        ! Compute velocities
        !
        u = rhou/rho
        v = rhov/rho
        w = rhow/rho


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
        ! Compute Kinetic Energy Jacobians
        !
        dke_drho  = -HALF*(u*u + v*v + w*w)
        dke_drhou = u
        dke_drhov = v
        dke_drhow = w


        !
        ! Compute Pressure Jacobians
        !
        dp_drho  = -(gam-ONE)*dke_drho
        dp_drhou = -(gam-ONE)*dke_drhou
        dp_drhov = -(gam-ONE)*dke_drhov
        dp_drhow = -(gam-ONE)*dke_drhow
        dp_drhoE =  (gam-ONE)   ! No negative sign

        
        !
        ! Compute Temperature Jacobians
        !
        const = ONE/287.15_rk
        dT_drho  = const*invrho*dp_drho  -  const*invrho*invrho*p
        dT_drhou = const*invrho*dp_drhou
        dT_drhov = const*invrho*dp_drhov
        dT_drhow = const*invrho*dp_drhow
        dT_drhoE = const*invrho*dp_drhoE



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
        ! Compute temperature
        !
        T = prop%fluid%compute_temperature(rho,rhou,rhov,rhow,rhoE)


        !
        ! Compute dynamic viscosity, second coefficient of viscosity
        !
        mu    = prop%fluid%compute_viscosity_dynamic(T)
        lamda = prop%fluid%compute_viscosity_second(mu,T)


        !
        ! Compute temperature gradient
        !
        dT_dx = dT_drho*drho_dx + dT_drhou*drhou_dx + dT_drhov*drhov_dx + dT_drhow*drhow_dx + dT_drhoE*drhoE_dx
        dT_dy = dT_drho*drho_dy + dT_drhou*drhou_dy + dT_drhov*drhov_dy + dT_drhow*drhow_dy + dT_drhoE*drhoE_dy
        dT_dz = dT_drho*drho_dz + dT_drhou*drhou_dz + dT_drhov*drhov_dz + dT_drhow*drhow_dz + dT_drhoE*drhoE_dz


        !
        ! Compute shear stress components
        !
        tau_xx = TWO*mu*du_dx  +  lamda*(du_dx + dv_dy + dw_dz)
        tau_yy = TWO*mu*dv_dy  +  lamda*(du_dx + dv_dy + dw_dz)
        tau_zz = TWO*mu*dw_dz  +  lamda*(du_dx + dv_dy + dw_dz)

        tau_xy = mu*(du_dy + dv_dx)
        tau_xz = mu*(du_dz + dw_dx)
        tau_yz = mu*(dw_dy + dv_dz)





        !=================================================
        ! Mass flux
        !=================================================


        !=================================================
        ! x-momentum flux
        !=================================================
        flux_x = -tau_xx
        flux_y = -tau_xy
        flux_z = -tau_xz

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhou, integrand)

        !=================================================
        ! y-momentum flux
        !=================================================
        flux_x = -tau_xy
        flux_y = -tau_yy
        flux_z = -tau_yz

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhov, integrand)

        !=================================================
        ! z-momentum flux
        !=================================================
        flux_x = -tau_xz
        flux_y = -tau_yz
        flux_z = -tau_zz

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhow, integrand)

        !=================================================
        ! Energy flux
        !=================================================
        flux_x = -(1003._rk*mu/0.8_rk)*dT_dx  -  (u*tau_xx + v*tau_xy + w*tau_xz)
        flux_y = -(1003._rk*mu/0.8_rk)*dT_dy  -  (u*tau_xy + v*tau_yy + w*tau_yz)
        flux_z = -(1003._rk*mu/0.8_rk)*dT_dz  -  (u*tau_xz + v*tau_yz + w*tau_zz)

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhoE, integrand)

    end subroutine compute
    !**********************************************************************************************























end module fluid_viscous_bc_operator
