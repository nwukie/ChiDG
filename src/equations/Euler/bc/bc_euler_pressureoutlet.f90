module bc_euler_pressureoutlet
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, HALF, ME

    use type_bc,                only: bc_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
    use DNAD_D
    
    use EULER_properties,       only: EULER_properties_t
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !----------------------------------------------------------------------------------------
    type, public, extends(bc_t) :: euler_pressureoutlet_t

    contains

        procedure   :: add_options  !< Add boundary condition options
        procedure   :: compute      !< boundary condition function implementation

    end type euler_pressureoutlet_t
    !****************************************************************************************




contains


    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_options(self)    
        class(euler_pressureoutlet_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('euler_pressureoutlet')


        !
        ! Add functions
        !
        call self%bcproperties%add('StaticPressure','Required')         ! add StaticPressure


        !
        ! Add parameters
        !


    end subroutine add_options
    !******************************************************************************************








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
        class(euler_pressureoutlet_t),  intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop


        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            rho_m,  rhou_m, rhov_m, rhow_m, rhoE_m, &
            flux_x, flux_y, flux_z, integrand,      &
            u_m,    v_m,    w_m,                    &
            H_bc,   rhoE_bc, gam_m


        real(rk)                                    :: time
        type(point_t),  allocatable, dimension(:)   :: coords
        real(rk),       allocatable, dimension(:)   ::  &
            p_bc, normx, normy, normz


        !
        ! Get equation indices
        !
        irho  = prop%get_eqn_index("rho")
        irhou = prop%get_eqn_index("rhou")
        irhov = prop%get_eqn_index("rhov")
        irhow = prop%get_eqn_index("rhow")
        irhoE = prop%get_eqn_index("rhoE")



        !
        ! Get back pressure from function.
        !
        coords = worker%coords()
        time   = worker%time()
        p_bc = self%bcproperties%compute("StaticPressure",time,coords)




        !
        ! Interpolate interior solution to face quadrature nodes
        !
        rho_m  = worker%interpolate(irho,  'value', ME)
        rhou_m = worker%interpolate(irhou, 'value', ME)
        rhov_m = worker%interpolate(irhov, 'value', ME)
        rhow_m = worker%interpolate(irhow, 'value', ME)
        rhoE_m = worker%interpolate(irhoE, 'value', ME)


        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)



        !
        ! Compute gamma
        !
        call prop%fluid%compute_gamma(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,gam_m)



        !
        ! Compute velocity components
        !
        u_m = rhou_m/rho_m
        v_m = rhov_m/rho_m
        w_m = rhow_m/rho_m



        !
        ! Compute boundary condition energy and enthalpy
        !
        rhoE_bc = p_bc/(gam_m - ONE) + (rho_m/TWO)*(u_m*u_m + v_m*v_m + w_m*w_m)
        H_bc = (rhoE_bc + p_bc)/rho_m




        !=================================================
        ! Mass flux
        !=================================================
        flux_x = (rho_m * u_m)
        flux_y = (rho_m * v_m)
        flux_z = (rho_m * w_m)

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irho, integrand)

        !=================================================
        ! x-momentum flux
        !=================================================
        flux_x = (rho_m * u_m * u_m) + p_bc
        flux_y = (rho_m * u_m * v_m)
        flux_z = (rho_m * u_m * w_m)

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhou, integrand)

        !=================================================
        ! y-momentum flux
        !=================================================
        flux_x = (rho_m * v_m * u_m)
        flux_y = (rho_m * v_m * v_m) + p_bc
        flux_z = (rho_m * v_m * w_m)

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhov, integrand)

        !=================================================
        ! z-momentum flux
        !=================================================
        flux_x = (rho_m * w_m * u_m)
        flux_y = (rho_m * w_m * v_m)
        flux_z = (rho_m * w_m * w_m) + p_bc

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhow, integrand)


        !=================================================
        ! Energy flux
        !=================================================
        flux_x = (rho_m * u_m * H_bc)
        flux_y = (rho_m * v_m * H_bc)
        flux_z = (rho_m * w_m * H_bc)

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhoE, integrand)


    end subroutine compute
    !**********************************************************************************************






end module bc_euler_pressureoutlet
