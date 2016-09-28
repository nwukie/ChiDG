module bc_state_pressureoutlet
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, HALF, ME

    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
    use DNAD_D
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !----------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: pressureoutlet_t

    contains

        procedure   :: init                 !< Set-up bc state with options/name etc.
        procedure   :: compute_bc_state     !< boundary condition function implementation

    end type pressureoutlet_t
    !****************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(pressureoutlet_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Pressure Outlet")


!        !
!        ! Set operator equations
!        !
!        call self%set_equation("Density"   )
!        call self%set_equation("X-Momentum")
!        call self%set_equation("Y-Momentum")
!        call self%set_equation("Z-Momentum")
!        call self%set_equation("Energy"    )


        !
        ! Add functions
        !
        call self%bcproperties%add('Static Pressure','Required')         ! add StaticPressure


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
    subroutine compute_bc_state(self,worker,prop)
        class(pressureoutlet_t),    intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Equation indices
        integer(ik) :: irho, irhou, irhov, irhow, irhoE


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::      &
            rho_m,  rhou_m,  rhov_m,  rhow_m,  rhoE_m,  &
            rho_bc, rhou_bc, rhov_bc, rhow_bc, rhoE_bc, &
            drho_dx_m, drhou_dx_m, drhov_dx_m, drhow_dx_m, drhoE_dx_m,  &
            drho_dy_m, drhou_dy_m, drhov_dy_m, drhow_dy_m, drhoE_dy_m,  &
            drho_dz_m, drhou_dz_m, drhov_dz_m, drhow_dz_m, drhoE_dz_m,  &
            flux_x, flux_y,  flux_z,  integrand,        &
            u_bc,   v_bc,    w_bc,                      &
            H_bc,   gam_m


        real(rk)                                    :: time
        type(point_t),  allocatable, dimension(:)   :: coords
        real(rk),       allocatable, dimension(:)   ::  &
            p_bc, normx, normy, normz


        !
        ! Get equation indices
        !
        irho  = prop%get_equation_index("Density"   )
        irhou = prop%get_equation_index("X-Momentum")
        irhov = prop%get_equation_index("Y-Momentum")
        irhow = prop%get_equation_index("Z-Momentum")
        irhoE = prop%get_equation_index("Energy"    )


        !
        ! Get back pressure from function.
        !
        coords = worker%coords()
        time   = worker%time()
        p_bc = self%bcproperties%compute("Static Pressure",time,coords)


        !
        ! Interpolate interior solution to face quadrature nodes
        !
        rho_m  = worker%interpolate(irho,  'value', ME)
        rhou_m = worker%interpolate(irhou, 'value', ME)
        rhov_m = worker%interpolate(irhov, 'value', ME)
        rhow_m = worker%interpolate(irhow, 'value', ME)
        rhoE_m = worker%interpolate(irhoE, 'value', ME)


        drho_dx_m  = worker%get_face_variable(irho,  'ddx', ME)
        drho_dy_m  = worker%get_face_variable(irho,  'ddy', ME)
        drho_dz_m  = worker%get_face_variable(irho,  'ddz', ME)

        drhou_dx_m = worker%get_face_variable(irhou, 'ddx', ME)
        drhou_dy_m = worker%get_face_variable(irhou, 'ddy', ME)
        drhou_dz_m = worker%get_face_variable(irhou, 'ddz', ME)

        drhov_dx_m = worker%get_face_variable(irhov, 'ddx', ME)
        drhov_dy_m = worker%get_face_variable(irhov, 'ddy', ME)
        drhov_dz_m = worker%get_face_variable(irhov, 'ddz', ME)

        drhow_dx_m = worker%get_face_variable(irhow, 'ddx', ME)
        drhow_dy_m = worker%get_face_variable(irhow, 'ddy', ME)
        drhow_dz_m = worker%get_face_variable(irhow, 'ddz', ME)
        
        drhoE_dx_m = worker%get_face_variable(irhoE, 'ddx', ME)
        drhoE_dy_m = worker%get_face_variable(irhoE, 'ddy', ME)
        drhoE_dz_m = worker%get_face_variable(irhoE, 'ddz', ME)










        !
        ! Compute gamma
        !
        gam_m = prop%fluid%compute_gamma(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m)


        !
        ! Extrapolate density and momentum
        !
        rho_bc  = rho_m
        rhou_bc = rhou_m
        rhov_bc = rhov_m
        rhow_bc = rhow_m


        !
        ! Compute velocities
        !
        u_bc = rhou_bc/rho_bc
        v_bc = rhov_bc/rho_bc
        w_bc = rhow_bc/rho_bc


        !
        ! Compute boundary condition energy and enthalpy
        !
        rhoE_bc = p_bc/(gam_m - ONE) + (rho_bc*HALF)*(u_bc*u_bc + v_bc*v_bc + w_bc*w_bc)


        !
        ! Store boundary condition state
        !
        call worker%store_bc_state(irho,  rho_bc, 'value')
        call worker%store_bc_state(irhou, rhou_bc,'value')
        call worker%store_bc_state(irhov, rhov_bc,'value')
        call worker%store_bc_state(irhow, rhow_bc,'value')
        call worker%store_bc_state(irhoE, rhoE_bc,'value')






        call worker%store_bc_state(irho, drho_dx_m, 'ddx')
        call worker%store_bc_state(irho, drho_dy_m, 'ddy')
        call worker%store_bc_state(irho, drho_dz_m, 'ddz')

        call worker%store_bc_state(irhou, drhou_dx_m, 'ddx')
        call worker%store_bc_state(irhou, drhou_dy_m, 'ddy')
        call worker%store_bc_state(irhou, drhou_dz_m, 'ddz')

        call worker%store_bc_state(irhov, drhov_dx_m, 'ddx')
        call worker%store_bc_state(irhov, drhov_dy_m, 'ddy')
        call worker%store_bc_state(irhov, drhov_dz_m, 'ddz')

        call worker%store_bc_state(irhow, drhow_dx_m, 'ddx')
        call worker%store_bc_state(irhow, drhow_dy_m, 'ddy')
        call worker%store_bc_state(irhow, drhow_dz_m, 'ddz')

        call worker%store_bc_state(irhoE, drhoE_dx_m, 'ddx')
        call worker%store_bc_state(irhoE, drhoE_dy_m, 'ddy')
        call worker%store_bc_state(irhoE, drhoE_dz_m, 'ddz')







    end subroutine compute_bc_state
    !**********************************************************************************************






end module bc_state_pressureoutlet
