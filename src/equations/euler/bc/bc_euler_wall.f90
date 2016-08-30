module bc_euler_wall
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO, ME

    use type_bc,                only: bc_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none
    


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_t) :: euler_wall_t

    contains
        procedure   :: add_options
        procedure   :: compute    !> bc implementation
    end type euler_wall_t
    !*******************************************************************************************




contains



    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_options(self)    
        class(euler_wall_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('euler_wall')


        !
        ! Add functions
        !


        !
        ! Add parameters
        !


    end subroutine add_options
    !******************************************************************************************















    !> Specialized compute routine for Euler Slip-Wall Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[in]      ielem   Index of the element being computed
    !!  @param[in]      iface   Index of the face being computed
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(euler_wall_t),            intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop


        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE

        real(rk)    :: gam_m

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            rho_m,  rhou_m, rhov_m, rhow_m, rhoE_m, p_m, integrand, flux_x, flux_y, flux_z,  &
            rhou_bc, rhov_bc, rhow_bc, rhoE_bc, u_bc, v_bc, w_bc, u_m, v_m, w_m, p_bc

        real(rk),   allocatable, dimension(:)   ::  &
            normx, normy, normz



        !
        ! Get equation indices
        !
        irho  = prop%get_eqn_index("rho")
        irhou = prop%get_eqn_index("rhou")
        irhov = prop%get_eqn_index("rhov")
        irhow = prop%get_eqn_index("rhow")
        irhoE = prop%get_eqn_index("rhoE")



        !
        ! Interpolate interior solution to quadrature nodes
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
        ! Compute interior pressure
        !
        call prop%fluid%compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,p_m)
        p_bc = p_m



        !
        ! Initialize arrays
        !
        flux_x = p_bc
        flux_y = p_bc
        flux_z = p_bc



        !
        ! Mass Flux
        !
        flux_x = ZERO
        flux_y = ZERO
        flux_z = ZERO
        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irho, integrand)


        !
        ! Add pressure flux to momentum equation
        !
        flux_x = p_bc
        flux_y = ZERO
        flux_z = ZERO
        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhou, integrand)



        !
        ! Add pressure flux to momentum equation
        !
        flux_x = ZERO
        flux_y = p_bc
        flux_z = ZERO

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhov, integrand)



        !
        ! Add pressure flux to momentum equation
        !
        flux_x = ZERO
        flux_y = ZERO
        flux_z = p_bc

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhow, integrand)


        !
        ! Energy Flux
        !
        flux_x = ZERO
        flux_y = ZERO
        flux_z = ZERO

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhoE, integrand)


    end subroutine compute
    !*****************************************************************************************************






end module bc_euler_wall
