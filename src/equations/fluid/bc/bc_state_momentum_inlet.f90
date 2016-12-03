module bc_state_momentum_inlet
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO

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
    !!  @date   2/8/2016
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: momentum_inlet_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type momentum_inlet_t
    !*******************************************************************************************




contains





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(momentum_inlet_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Momentum Inlet")
        call self%set_family("Inlet")


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
        call self%bcproperties%add('Density',   'Required')
        call self%bcproperties%add('X-Velocity','Required')
        call self%bcproperties%add('Y-Velocity','Required')
        call self%bcproperties%add('Z-Velocity','Required')



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
        class(momentum_inlet_t),    intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        class(properties_t),    intent(inout)   :: prop



        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::              &
            rho_m,  rhou_m,  rhov_m,  rhow_m,  rhoE_m,  p_m,    &
            rho_bc, rhou_bc, rhov_bc, rhow_bc, rhoE_bc, p_bc,   &
            drho_dx_m, drhou_dx_m, drhov_dx_m, drhow_dx_m, drhoE_dx_m,  &
            drho_dy_m, drhou_dy_m, drhov_dy_m, drhow_dy_m, drhoE_dy_m,  &
            drho_dz_m, drhou_dz_m, drhov_dz_m, drhow_dz_m, drhoE_dz_m,  &
            flux_x, flux_y, flux_z, integrand,                  &
            u_m,    v_m,    w_m,                                &
            u_bc,   v_bc,   w_bc


        real(rk)                                    :: time, gam_m
        type(point_t),  allocatable, dimension(:)   :: coords



        !
        ! Get equation indices
        !
        irho  = prop%get_primary_field_index("Density"   )
        irhou = prop%get_primary_field_index("X-Momentum")
        irhov = prop%get_primary_field_index("Y-Momentum")
        irhow = prop%get_primary_field_index("Z-Momentum")
        irhoE = prop%get_primary_field_index("Energy"    )




        !
        ! Interpolate interior solution to quadrature nodes
        !
        rho_m  = worker%get_primary_field_face("Density"   , 'value', 'face interior')
        rhou_m = worker%get_primary_field_face("X-Momentum", 'value', 'face interior')
        rhov_m = worker%get_primary_field_face("Y-Momentum", 'value', 'face interior')
        rhow_m = worker%get_primary_field_face("Z-Momentum", 'value', 'face interior')
        rhoE_m = worker%get_primary_field_face("Energy"    , 'value', 'face interior')

        drho_dx_m  = worker%get_primary_field_face("Density"   , 'ddx', 'face interior')
        drho_dy_m  = worker%get_primary_field_face("Density"   , 'ddy', 'face interior')
        drho_dz_m  = worker%get_primary_field_face("Density"   , 'ddz', 'face interior')

        drhou_dx_m = worker%get_primary_field_face("X-Momentum", 'ddx', 'face interior')
        drhou_dy_m = worker%get_primary_field_face("X-Momentum", 'ddy', 'face interior')
        drhou_dz_m = worker%get_primary_field_face("X-Momentum", 'ddz', 'face interior')

        drhov_dx_m = worker%get_primary_field_face("Y-Momentum", 'ddx', 'face interior')
        drhov_dy_m = worker%get_primary_field_face("Y-Momentum", 'ddy', 'face interior')
        drhov_dz_m = worker%get_primary_field_face("Y-Momentum", 'ddz', 'face interior')

        drhow_dx_m = worker%get_primary_field_face("Z-Momentum", 'ddx', 'face interior')
        drhow_dy_m = worker%get_primary_field_face("Z-Momentum", 'ddy', 'face interior')
        drhow_dz_m = worker%get_primary_field_face("Z-Momentum", 'ddz', 'face interior')
        
        drhoE_dx_m = worker%get_primary_field_face("Energy"    , 'ddx', 'face interior')
        drhoE_dy_m = worker%get_primary_field_face("Energy"    , 'ddy', 'face interior')
        drhoE_dz_m = worker%get_primary_field_face("Energy"    , 'ddz', 'face interior')


        !
        ! Initialize variables
        !
        rho_bc = rho_m
        u_bc   = rho_m
        v_bc   = rho_m
        w_bc   = rho_m


        !
        ! Get boundary condition Total Temperature, Total Pressure, and normal vector
        !
        coords = worker%coords()
        time   = worker%time()
        rho_bc = self%bcproperties%compute("Density",     time, coords)
        u_bc   = self%bcproperties%compute("X-Velocity",  time, coords)
        v_bc   = self%bcproperties%compute("Y-Velocity",  time, coords)
        w_bc   = self%bcproperties%compute("Z-Velocity",  time, coords)





        !
        ! Compute bc momentum
        !
        rhou_bc = rho_bc * u_bc
        rhov_bc = rho_bc * v_bc
        rhow_bc = rho_bc * w_bc

        !
        ! Compute interior pressure
        !
        !p_m = prop%fluid%compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m)
        !gam_m = prop%fluid%compute_gamma(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m)
        p_m = worker%get_model_field_face('Pressure', 'value', 'face interior')
        gam_m = 1.4_rk



        !
        ! Compute bc energy
        !
        rhoE_bc = p_m/(gam_m - ONE) + (rho_bc/TWO)*( (u_bc*u_bc) + (v_bc*v_bc) + (w_bc*w_bc) )



        !
        ! Store computed boundary state
        !
        call worker%store_bc_state("Density"   ,rho_bc, 'value')
        call worker%store_bc_state("X-Momentum",rhou_bc,'value')
        call worker%store_bc_state("Y-Momentum",rhov_bc,'value')
        call worker%store_bc_state("Z-Momentum",rhow_bc,'value')
        call worker%store_bc_state("Energy"    ,rhoE_bc,'value')




        call worker%store_bc_state("Density"   , drho_dx_m,  'ddx')
        call worker%store_bc_state("Density"   , drho_dy_m,  'ddy')
        call worker%store_bc_state("Density"   , drho_dz_m,  'ddz')

        call worker%store_bc_state("X-Momentum", drhou_dx_m, 'ddx')
        call worker%store_bc_state("X-Momentum", drhou_dy_m, 'ddy')
        call worker%store_bc_state("X-Momentum", drhou_dz_m, 'ddz')
                                                
        call worker%store_bc_state("Y-Momentum", drhov_dx_m, 'ddx')
        call worker%store_bc_state("Y-Momentum", drhov_dy_m, 'ddy')
        call worker%store_bc_state("Y-Momentum", drhov_dz_m, 'ddz')
                                                
        call worker%store_bc_state("Z-Momentum", drhow_dx_m, 'ddx')
        call worker%store_bc_state("Z-Momentum", drhow_dy_m, 'ddy')
        call worker%store_bc_state("Z-Momentum", drhow_dz_m, 'ddz')
                                                
        call worker%store_bc_state("Energy"    , drhoE_dx_m, 'ddx')
        call worker%store_bc_state("Energy"    , drhoE_dy_m, 'ddy')
        call worker%store_bc_state("Energy"    , drhoE_dz_m, 'ddz')










    end subroutine compute_bc_state
    !******************************************************************************************









end module bc_state_momentum_inlet
