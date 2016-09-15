module bc_state_totalinlet
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, ME

    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
    use DNAD_D
    
!    use EULER_properties,       only: EULER_properties_t
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: totalinlet_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type totalinlet_t
    !*******************************************************************************************




contains





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(totalinlet_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Total Inlet")


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

    end subroutine init
    !********************************************************************************










    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   9/8/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop)
        class(totalinlet_t),    intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        class(properties_t),    intent(inout)   :: prop



        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::              &
            rho_m,  rhou_m,  rhov_m,  rhow_m,  rhoE_m,  p_m,    &
            rho_bc, rhou_bc, rhov_bc, rhow_bc, rhoE_bc, p_bc,   &
            flux_x, flux_y, flux_z, integrand,                  &
            u_m,    v_m,    w_m,                                &
            u_bc,   v_bc,   w_bc,                               &
            T_bc,   vmag2_m, vmag, H_bc


        real(rk)                                    :: gam_m, cp_m, M, time
        type(point_t),  allocatable, dimension(:)   :: coords
        real(rk),       allocatable, dimension(:)   ::  &
            TT, PT, nx, ny, nz, nr, nt, theta,          &
            normal_direction, x, y, r, normx, normy, normz



        !
        ! Get equation indices
        !
        irho  = prop%get_equation_index("Density"   )
        irhou = prop%get_equation_index("X-Momentum")
        irhov = prop%get_equation_index("Y-Momentum")
        irhow = prop%get_equation_index("Z-Momentum")
        irhoE = prop%get_equation_index("Energy"    )


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

!        select type(prop)
!            type is (EULER_properties_t)
!                cp_m  = (prop%R)*(gam_m/(gam_m-ONE))
!        end select
         cp_m  = 287.15_rk*(gam_m/(gam_m-ONE))


        T_bc = TT - (vmag2_m)/(TWO*cp_m)
        p_bc = PT*((T_bc/TT)**(gam_m/(gam_m-ONE)))


        !
        ! Compute boundary condition density from ideal gas law
        !
!        select type(prop)
!            type is (EULER_properties_t)
!                rho_bc = p_bc/(T_bc*prop%R)
!        end select
         rho_bc = p_bc/(T_bc*287.15_rk)


        !
        ! Compute bc momentum
        !
        rhou_bc = rho_bc * u_bc
        rhov_bc = rho_bc * v_bc
        rhow_bc = rho_bc * w_bc


        !
        ! Compute bc energy
        !
        rhoE_bc = p_bc/(gam_m - ONE) + (rho_bc/TWO)*( (u_bc*u_bc) + (v_bc*v_bc) + (w_bc*w_bc) )



        !
        ! Store computed boundary state
        !
        call worker%store_bc_state(irho, rho_bc, 'value')
        call worker%store_bc_state(irhou,rhou_bc,'value')
        call worker%store_bc_state(irhov,rhov_bc,'value')
        call worker%store_bc_state(irhow,rhow_bc,'value')
        call worker%store_bc_state(irhoE,rhoE_bc,'value')


    end subroutine compute_bc_state
    !******************************************************************************************









end module bc_state_totalinlet
