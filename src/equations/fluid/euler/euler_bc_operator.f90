module euler_bc_operator
    use mod_kinds,          only: ik, rk
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
    type, public, extends(operator_t)   :: euler_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type euler_bc_operator_t
    !*************************************************************************************






contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Euler BC Flux")

        !
        ! Set operator type
        !
        call self%set_operator_type("BC Advective Flux")

        !
        ! Set operator equations
        !
        call self%add_primary_field("Density"   )
        call self%add_primary_field("X-Momentum")
        call self%add_primary_field("Y-Momentum")
        call self%add_primary_field("Z-Momentum")
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
        class(euler_bc_operator_t), intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::          &
            rho_bc,  rhou_bc, rhov_bc, rhow_bc, rhoE_bc,    &
            u_bc,    v_bc,    w_bc,                         &
            H_bc,    gam_bc,  p_bc,                         &
            flux_x,  flux_y,  flux_z,  integrand


        real(rk),   allocatable, dimension(:)   ::          &
            normx, normy, normz


        !
        ! Get equation indices
        !
        irho  = prop%get_primary_field_index("Density"   )
        irhou = prop%get_primary_field_index("X-Momentum")
        irhov = prop%get_primary_field_index("Y-Momentum")
        irhow = prop%get_primary_field_index("Z-Momentum")
        irhoE = prop%get_primary_field_index("Energy"    )



        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        rho_bc  = worker%get_primary_field_face("Density"   ,irho,  'value', 'boundary')
        rhou_bc = worker%get_primary_field_face("X-Momentum",irhou, 'value', 'boundary')
        rhov_bc = worker%get_primary_field_face("Y-Momentum",irhov, 'value', 'boundary')
        rhow_bc = worker%get_primary_field_face("Z-Momentum",irhow, 'value', 'boundary')
        rhoE_bc = worker%get_primary_field_face("Energy"    ,irhoE, 'value', 'boundary')


        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)



        !
        ! Compute gamma
        !
        gam_bc = prop%fluid%compute_gamma(rho_bc,rhou_bc,rhov_bc,rhow_bc,rhoE_bc)
        p_bc   = prop%fluid%compute_pressure(rho_bc,rhou_bc,rhov_bc,rhow_bc,rhoE_bc)



        !
        ! Compute velocity components
        !
        u_bc = rhou_bc/rho_bc
        v_bc = rhov_bc/rho_bc
        w_bc = rhow_bc/rho_bc



        !
        ! Compute boundary condition energy and enthalpy
        !
        H_bc = (rhoE_bc + p_bc)/rho_bc




        !=================================================
        ! Mass flux
        !=================================================
        flux_x = (rho_bc * u_bc)
        flux_y = (rho_bc * v_bc)
        flux_z = (rho_bc * w_bc)

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary('Density',irho, integrand)

        !=================================================
        ! x-momentum flux
        !=================================================
        flux_x = (rho_bc * u_bc * u_bc) + p_bc
        flux_y = (rho_bc * u_bc * v_bc)
        flux_z = (rho_bc * u_bc * w_bc)

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary('X-Momentum',irhou, integrand)

        !=================================================
        ! y-momentum flux
        !=================================================
        flux_x = (rho_bc * v_bc * u_bc)
        flux_y = (rho_bc * v_bc * v_bc) + p_bc
        flux_z = (rho_bc * v_bc * w_bc)

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary('Y-Momentum',irhov, integrand)

        !=================================================
        ! z-momentum flux
        !=================================================
        flux_x = (rho_bc * w_bc * u_bc)
        flux_y = (rho_bc * w_bc * v_bc)
        flux_z = (rho_bc * w_bc * w_bc) + p_bc

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary('Z-Momentum',irhow, integrand)

        !=================================================
        ! Energy flux
        !=================================================
        flux_x = (rho_bc * u_bc * H_bc)
        flux_y = (rho_bc * v_bc * H_bc)
        flux_z = (rho_bc * w_bc * H_bc)

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary('Energy',irhoE, integrand)

    end subroutine compute
    !**********************************************************************************************























end module euler_bc_operator
