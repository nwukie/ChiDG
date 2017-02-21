module fluid_viscous_bc_operator
    use mod_kinds,          only: ik, rk
    use mod_constants,      only: HALF, ONE, TWO
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
        call self%add_primary_field("Density"   )
        call self%add_primary_field("Momentum-1")
        call self%add_primary_field("Momentum-2")
        call self%add_primary_field("Momentum-3")
        call self%add_primary_field("Energy"    )

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


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::                &
            rho, rhou, rhov, rhow, rhoE, p, u, v, w, invrho,    &
            drho_dx, drhou_dx, drhov_dx, drhow_dx, drhoE_dx,    &
            drho_dy, drhou_dy, drhov_dy, drhow_dy, drhoE_dy,    &
            drho_dz, drhou_dz, drhov_dz, drhow_dz, drhoE_dz,    &
            dT_drho, dT_drhou, dT_drhov, dT_drhow, dT_drhoE,    &
            dp_drho, dp_drhou, dp_drhov, dp_drhow, dp_drhoE,    &
            dke_drho, dke_drhou, dke_drhov, dke_drhow,          &
            dT_dx,   dT_dy,    dT_dz,                           &
            k,       k_l,      k_t,                             &
            tau_xx, tau_yy, tau_zz, tau_xy, tau_xz, tau_yz,     &
            flux_x, flux_y, flux_z, integrand

        real(rk),   allocatable, dimension(:)   ::          &
            normx, normy, normz

        real(rk) :: const, gam


        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        rho  = worker%get_primary_field_face("Density"   ,'value', 'boundary')
        rhou = worker%get_primary_field_face("Momentum-1",'value', 'boundary')
        rhov = worker%get_primary_field_face("Momentum-2",'value', 'boundary')
        rhow = worker%get_primary_field_face("Momentum-3",'value', 'boundary')
        rhoE = worker%get_primary_field_face("Energy"    ,'value', 'boundary')


        !
        ! Interpolate solution gradients to quadrature nodes
        !
        drho_dx  = worker%get_primary_field_face("Density"   ,'grad1+lift', 'boundary')
        drho_dy  = worker%get_primary_field_face("Density"   ,'grad2+lift', 'boundary')
        drho_dz  = worker%get_primary_field_face("Density"   ,'grad3+lift', 'boundary')

        drhou_dx = worker%get_primary_field_face("Momentum-1",'grad1+lift', 'boundary')
        drhou_dy = worker%get_primary_field_face("Momentum-1",'grad2+lift', 'boundary')
        drhou_dz = worker%get_primary_field_face("Momentum-1",'grad3+lift', 'boundary')

        drhov_dx = worker%get_primary_field_face("Momentum-2",'grad1+lift', 'boundary')
        drhov_dy = worker%get_primary_field_face("Momentum-2",'grad2+lift', 'boundary')
        drhov_dz = worker%get_primary_field_face("Momentum-2",'grad3+lift', 'boundary')

        drhow_dx = worker%get_primary_field_face("Momentum-3",'grad1+lift', 'boundary')
        drhow_dy = worker%get_primary_field_face("Momentum-3",'grad2+lift', 'boundary')
        drhow_dz = worker%get_primary_field_face("Momentum-3",'grad3+lift', 'boundary')

        drhoE_dx = worker%get_primary_field_face("Energy"    ,'grad1+lift', 'boundary')
        drhoE_dy = worker%get_primary_field_face("Energy"    ,'grad2+lift', 'boundary')
        drhoE_dz = worker%get_primary_field_face("Energy"    ,'grad3+lift', 'boundary')





        !
        ! Get normal vector
        !
        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)


        !
        ! Get Model fields:
        !   Pressure
        !   Thermal Conductivity
        !
        p       = worker%get_model_field_face('Pressure',                       'value', 'boundary')
        k_l     = worker%get_model_field_face('Laminar Thermal Conductivity',   'value', 'boundary')
        k_t     = worker%get_model_field_face('Turbulent Thermal Conductivity', 'value', 'boundary')
        gam = 1.4_rk



        !
        ! Compute effective conductivity. Laminar + Turbulent.
        !
        k = k_l + k_t





        !
        ! Compute velocities
        !
        invrho = ONE/rho
        u = rhou*invrho
        v = rhov*invrho
        w = rhow*invrho



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
        dp_drhoE =  dp_drhow    ! Initialize derivatives
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
        ! Compute temperature gradient
        !
        dT_dx = dT_drho*drho_dx + dT_drhou*drhou_dx + dT_drhov*drhov_dx + dT_drhow*drhow_dx + dT_drhoE*drhoE_dx
        dT_dy = dT_drho*drho_dy + dT_drhou*drhou_dy + dT_drhov*drhov_dy + dT_drhow*drhow_dy + dT_drhoE*drhoE_dy
        dT_dz = dT_drho*drho_dz + dT_drhou*drhou_dz + dT_drhov*drhov_dz + dT_drhow*drhow_dz + dT_drhoE*drhoE_dz


        !
        ! get shear stress components
        !
        tau_xx = worker%get_model_field_face('Shear-11', 'value', 'boundary')
        tau_yy = worker%get_model_field_face('Shear-22', 'value', 'boundary')
        tau_zz = worker%get_model_field_face('Shear-33', 'value', 'boundary')

        tau_xy = worker%get_model_field_face('Shear-12', 'value', 'boundary')
        tau_xz = worker%get_model_field_face('Shear-13', 'value', 'boundary')
        tau_yz = worker%get_model_field_face('Shear-23', 'value', 'boundary')

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

        call worker%integrate_boundary('Momentum-1',integrand)

        !=================================================
        ! y-momentum flux
        !=================================================
        flux_x = -tau_xy
        flux_y = -tau_yy
        flux_z = -tau_yz

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary('Momentum-2',integrand)

        !=================================================
        ! z-momentum flux
        !=================================================
        flux_x = -tau_xz
        flux_y = -tau_yz
        flux_z = -tau_zz

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary('Momentum-3',integrand)

        !=================================================
        ! Energy flux
        !=================================================
        flux_x = -k*dT_dx  -  (u*tau_xx + v*tau_xy + w*tau_xz)
        flux_y = -k*dT_dy  -  (u*tau_xy + v*tau_yy + w*tau_yz)
        flux_z = -k*dT_dz  -  (u*tau_xz + v*tau_yz + w*tau_zz)

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary('Energy',integrand)

    end subroutine compute
    !**********************************************************************************************























end module fluid_viscous_bc_operator
