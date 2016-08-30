module bc_euler_totalinlet
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, HALF, ZERO, ME

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
    !!  @date   2/8/2016
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_t) :: euler_totalinlet_t

    contains

        procedure   :: add_options  !< Add boundary condition options
        procedure   :: compute      !< bc implementation

    end type euler_totalinlet_t
    !*******************************************************************************************




contains




    !>  Add options for total pressure/temperature boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine add_options(self)
        class(euler_totalinlet_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('euler_totalinlet')


        !
        ! Add functions
        !
        call self%bcproperties%add('TotalPressure',   'Required')
        call self%bcproperties%add('TotalTemperature','Required')

        call self%bcproperties%add('normal_direction','Required')
        call self%bcproperties%add('nx',              'Required')
        call self%bcproperties%add('ny',              'Required')
        call self%bcproperties%add('nz',              'Required')
        call self%bcproperties%add('nr',              'Required')
        call self%bcproperties%add('nt',              'Required')

        !
        ! Set default angle
        !
        call self%set_fcn_option('nx', 'val', 1._rk)
        call self%set_fcn_option('ny', 'val', 0._rk)
        call self%set_fcn_option('nz', 'val', 0._rk)


    end subroutine add_options
    !********************************************************************************************







    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/6/2016
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[in]      ielem   Index of the element being computed
    !!  @param[in]      iface   Index of the face being computed
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(euler_totalinlet_t),      intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop


        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            rho_m,  rhou_m, rhov_m, rhow_m, rhoE_m, p_m,        &
            flux_x, flux_y, flux_z, integrand,                  &
            u_m,    v_m,    w_m,                                &
            u_bc,   v_bc,   w_bc,                               &
            T_bc,   p_bc,   rho_bc, rhoE_bc,                    &
            vmag2_m, vmag, H_bc


        real(rk)                                    :: gam_m, cp_m, M, time
        type(point_t),  allocatable, dimension(:)   :: coords
        real(rk),       allocatable, dimension(:)   ::  &
            TT, PT, nx, ny, nz, nr, nt, theta,          &
            normal_direction, x, y, r, normx, normy, normz




        !
        ! Get equation indices
        !
        irho  = prop%get_eqn_index("rho")
        irhou = prop%get_eqn_index("rhou")
        irhov = prop%get_eqn_index("rhov")
        irhow = prop%get_eqn_index("rhow")
        irhoE = prop%get_eqn_index("rhoE")


        !
        ! Get boundary condition Total Temperature, Total Pressure, and normal vector
        !
        coords = worker%coords()
        time   = worker%time()
        PT = self%bcproperties%compute("TotalPressure",     time, coords)
        TT = self%bcproperties%compute("TotalTemperature",  time, coords)


        normal_direction = self%bcproperties%compute("normal_direction",    time, coords)
        nx = self%bcproperties%compute("nx",                                time, coords)
        ny = self%bcproperties%compute("ny",                                time, coords)
        nz = self%bcproperties%compute("nz",                                time, coords)



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
        ! Compute velocity components
        !
        u_m = rhou_m/rho_m
        v_m = rhov_m/rho_m
        w_m = rhow_m/rho_m


        !
        ! Compute velocity magnitude squared from interior state
        !
        vmag2_m = (u_m*u_m) + (v_m*v_m) + (w_m*w_m)
        vmag = sqrt(vmag2_m)


        !
        ! Compute boundary condition velocity components from imposed direction
        !
        u_bc = vmag*nx
        v_bc = vmag*ny
        w_bc = vmag*nz








        !
        ! Compute boundary condition temperature and pressure
        !
        !& HARDCODED GAMMA. HARDCODED CP
        gam_m = 1.4_rk

        select type(prop)
            type is (EULER_properties_t)
                cp_m  = (prop%R)*(gam_m/(gam_m-ONE))
        end select

        T_bc = TT - (vmag2_m)/(TWO*cp_m)
        p_bc = PT*((T_bc/TT)**(gam_m/(gam_m-ONE)))


        !
        ! Compute boundary condition density from ideal gas law
        !
        select type(prop)
            type is (EULER_properties_t)
                rho_bc = p_bc/(T_bc*prop%R)
        end select




        !
        ! Compute energy and enthalpy
        !
        rhoE_bc = p_bc/(gam_m - ONE) + (rho_bc/TWO)*( (u_bc*u_bc) + (v_bc*v_bc) + (w_bc*w_bc) )
        H_bc    = (rhoE_bc + p_bc)/rho_bc



        !=================================================
        ! Mass flux
        !=================================================
        flux_x = (rho_bc * u_bc)
        flux_y = (rho_bc * v_bc)
        flux_z = (rho_bc * w_bc)

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irho, integrand)

        !=================================================
        ! x-momentum flux
        !=================================================
        flux_x = (rho_bc * u_bc * u_bc) + p_bc
        flux_y = (rho_bc * u_bc * v_bc)
        flux_z = (rho_bc * u_bc * w_bc)

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhou, integrand)

        !=================================================
        ! y-momentum flux
        !=================================================
        flux_x = (rho_bc * v_bc * u_bc)
        flux_y = (rho_bc * v_bc * v_bc) + p_bc
        flux_z = (rho_bc * v_bc * w_bc)

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhov, integrand)

        !=================================================
        ! z-momentum flux
        !=================================================
        flux_x = (rho_bc * w_bc * u_bc)
        flux_y = (rho_bc * w_bc * v_bc)
        flux_z = (rho_bc * w_bc * w_bc) + p_bc

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhow, integrand)


        !=================================================
        ! Energy flux
        !=================================================
        flux_x = (rho_bc * u_bc * H_bc)
        flux_y = (rho_bc * v_bc * H_bc)
        flux_z = (rho_bc * w_bc * H_bc)

        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        call worker%integrate_boundary(irhoE, integrand)


    end subroutine compute
    !***********************************************************************************************






end module bc_euler_totalinlet
