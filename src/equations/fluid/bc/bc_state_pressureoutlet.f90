module bc_state_pressureoutlet
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, HALF

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
        ! Set name, family
        !
        call self%set_name("Pressure Outlet")
        call self%set_family("Outlet")


!        !
!        ! Set operator equations
!        !
!        call self%set_equation("Density"   )
!        call self%set_equation("Momentum-1")
!        call self%set_equation("Momentum-2")
!        call self%set_equation("Momentum-3")
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
            H_bc


        real(rk)                                    :: time, gam_m
        type(point_t),  allocatable, dimension(:)   :: coords
        real(rk),       allocatable, dimension(:)   ::  &
            p_bc, normx, normy, normz


        !
        ! Get equation indices
        !
        irho  = prop%get_primary_field_index("Density"   )
        irhou = prop%get_primary_field_index("Momentum-1")
        irhov = prop%get_primary_field_index("Momentum-2")
        irhow = prop%get_primary_field_index("Momentum-3")
        irhoE = prop%get_primary_field_index("Energy"    )


        !
        ! Get back pressure from function.
        !
        coords = worker%coords()
        time   = worker%time()
        p_bc = self%bcproperties%compute("Static Pressure",time,coords)


        !
        ! Interpolate interior solution to face quadrature nodes
        !
        !rho_m  = worker%interpolate(irho,  'value', 'face interior')
        !rhou_m = worker%interpolate(irhou, 'value', 'face interior')
        !rhov_m = worker%interpolate(irhov, 'value', 'face interior')
        !rhow_m = worker%interpolate(irhow, 'value', 'face interior')
        !rhoE_m = worker%interpolate(irhoE, 'value', 'face interior')
        rho_m  = worker%get_primary_field_face("Density"   ,'value', 'face interior')
        rhou_m = worker%get_primary_field_face("Momentum-1",'value', 'face interior')
        rhov_m = worker%get_primary_field_face("Momentum-2",'value', 'face interior')
        rhow_m = worker%get_primary_field_face("Momentum-3",'value', 'face interior')
        rhoE_m = worker%get_primary_field_face("Energy"    ,'value', 'face interior')


        drho_dx_m  = worker%get_primary_field_face("Density"   ,'grad1', 'face interior')
        drho_dy_m  = worker%get_primary_field_face("Density"   ,'grad2', 'face interior')
        drho_dz_m  = worker%get_primary_field_face("Density"   ,'grad3', 'face interior')

        drhou_dx_m = worker%get_primary_field_face("Momentum-1",'grad1', 'face interior')
        drhou_dy_m = worker%get_primary_field_face("Momentum-1",'grad2', 'face interior')
        drhou_dz_m = worker%get_primary_field_face("Momentum-1",'grad3', 'face interior')

        drhov_dx_m = worker%get_primary_field_face("Momentum-2",'grad1', 'face interior')
        drhov_dy_m = worker%get_primary_field_face("Momentum-2",'grad2', 'face interior')
        drhov_dz_m = worker%get_primary_field_face("Momentum-2",'grad3', 'face interior')

        drhow_dx_m = worker%get_primary_field_face("Momentum-3",'grad1', 'face interior')
        drhow_dy_m = worker%get_primary_field_face("Momentum-3",'grad2', 'face interior')
        drhow_dz_m = worker%get_primary_field_face("Momentum-3",'grad3', 'face interior')
        
        drhoE_dx_m = worker%get_primary_field_face("Energy"    ,'grad1', 'face interior')
        drhoE_dy_m = worker%get_primary_field_face("Energy"    ,'grad2', 'face interior')
        drhoE_dz_m = worker%get_primary_field_face("Energy"    ,'grad3', 'face interior')










        !
        ! Compute gamma
        !
        !gam_m = prop%fluid%compute_gamma(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m)
        gam_m = 1.4_rk


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
        call worker%store_bc_state("Density"   ,rho_bc, 'value')
        call worker%store_bc_state("Momentum-1",rhou_bc,'value')
        call worker%store_bc_state("Momentum-2",rhov_bc,'value')
        call worker%store_bc_state("Momentum-3",rhow_bc,'value')
        call worker%store_bc_state("Energy"    ,rhoE_bc,'value')





        call worker%store_bc_state("Density"   , drho_dx_m, 'grad1')
        call worker%store_bc_state("Density"   , drho_dy_m, 'grad2')
        call worker%store_bc_state("Density"   , drho_dz_m, 'grad3')
                                                
        call worker%store_bc_state("Momentum-1", drhou_dx_m, 'grad1')
        call worker%store_bc_state("Momentum-1", drhou_dy_m, 'grad2')
        call worker%store_bc_state("Momentum-1", drhou_dz_m, 'grad3')
                                                
        call worker%store_bc_state("Momentum-2", drhov_dx_m, 'grad1')
        call worker%store_bc_state("Momentum-2", drhov_dy_m, 'grad2')
        call worker%store_bc_state("Momentum-2", drhov_dz_m, 'grad3')
                                                
        call worker%store_bc_state("Momentum-3", drhow_dx_m, 'grad1')
        call worker%store_bc_state("Momentum-3", drhow_dy_m, 'grad2')
        call worker%store_bc_state("Momentum-3", drhow_dz_m, 'grad3')
                                                
        call worker%store_bc_state("Energy"    , drhoE_dx_m, 'grad1')
        call worker%store_bc_state("Energy"    , drhoE_dy_m, 'grad2')
        call worker%store_bc_state("Energy"    , drhoE_dz_m, 'grad3')







    end subroutine compute_bc_state
    !**********************************************************************************************






end module bc_state_pressureoutlet
