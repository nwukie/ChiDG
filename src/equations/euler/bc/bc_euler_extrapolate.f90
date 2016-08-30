module bc_euler_extrapolate
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, HALF, ZERO, ME

    use type_bc,                only: bc_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    
    use EULER_properties,       only: EULER_properties_t
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_t) :: euler_extrapolate_t

    contains

        procedure :: compute    !> bc implementation

    end type euler_extrapolate_t
    !-------------------------------------------------------------------------------------------




contains

    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[in]      ielem   Index of the element being computed
    !!  @param[in]      iface   Index of the face being computed
    !!  @param[in]      iblk    Index of the linearization block being computed
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(euler_extrapolate_t),     intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop


        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::          &
            rho_m,  rhou_m, rhov_m, rhow_m, rhoE_m, p_m,    &
            H_m,    u_m,    v_m,    w_m,                    &
            flux_x, flux_y, flux_z, integrand

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


        call prop%fluid%compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,p_m)

        u_m = rhou_m/rho_m
        v_m = rhov_m/rho_m
        w_m = rhow_m/rho_m

        H_m = (rhoE_m + p_m)/rho_m

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
        flux_x = (rho_m * u_m * u_m) + p_m
        flux_y = (rho_m * u_m * v_m)
        flux_z = (rho_m * u_m * w_m)

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhou, integrand)

        !=================================================
        ! y-momentum flux
        !=================================================
        flux_x = (rho_m * v_m * u_m)
        flux_y = (rho_m * v_m * v_m) + p_m
        flux_z = (rho_m * v_m * w_m)

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhov, integrand)

        !=================================================
        ! z-momentum flux
        !=================================================
        flux_x = (rho_m * w_m * u_m)
        flux_y = (rho_m * w_m * v_m)
        flux_z = (rho_m * w_m * w_m) + p_m

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhow, integrand)


        !=================================================
        ! Energy flux
        !=================================================
        flux_x = (rho_m * u_m * H_m)
        flux_y = (rho_m * v_m * H_m)
        flux_z = (rho_m * w_m * H_m)

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhoE, integrand)



    end subroutine compute
    !*************************************************************************************************






end module bc_euler_extrapolate
