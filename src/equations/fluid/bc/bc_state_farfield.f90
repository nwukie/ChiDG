module bc_state_farfield
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO, ONE, RKTOL, FOUR

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
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_state_t) :: farfield_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type farfield_t
    !*******************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(farfield_t),   intent(inout) :: self
        

        !
        ! Set operator name, family
        !
        call self%set_name("Farfield")
        call self%set_family("Farfield")



        call self%bcproperties%add('Density',   'Required')
        call self%bcproperties%add('Pressure',  'Required')
        call self%bcproperties%add('X-Velocity','Required')
        call self%bcproperties%add('Y-Velocity','Required')
        call self%bcproperties%add('Z-Velocity','Required')




    end subroutine init
    !********************************************************************************














    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop)
        class(farfield_t),      intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        class(properties_t),    intent(inout)   :: prop

        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::      &
            rho_m,  rhou_m,  rhov_m,  rhow_m,  rhoE_m,  &
            rho_bc, rhou_bc, rhov_bc, rhow_bc, rhoE_bc, p_bc, &
            drho_dx_m, drhou_dx_m, drhov_dx_m, drhow_dx_m, drhoE_dx_m,  &
            drho_dy_m, drhou_dy_m, drhov_dy_m, drhow_dy_m, drhoE_dy_m,  &
            drho_dz_m, drhou_dz_m, drhov_dz_m, drhow_dz_m, drhoE_dz_m,  &
            normal_momentum, normal_velocity, R_inf, R_extrapolated, &
            u_bc_norm, v_bc_norm, w_bc_norm, u_bc_tang, v_bc_tang, w_bc_tang, entropy_bc, &
            c_bc, c_m, p_m, T_m, u_m, v_m, w_m

        real(rk)    :: time

        real(rk), allocatable, dimension(:) ::              &
            unormx, unormy, unormz,                         &
            rho_input, p_input, u_input, v_input, w_input, T_input, c_input

        type(point_t),  allocatable, dimension(:)   :: coords

        logical, allocatable, dimension(:)  :: inflow, outflow



        !
        ! Get equation indices
        !
        irho  = prop%get_primary_field_index("Density"   )
        irhou = prop%get_primary_field_index("X-Momentum")
        irhov = prop%get_primary_field_index("Y-Momentum")
        irhow = prop%get_primary_field_index("Z-Momentum")
        irhoE = prop%get_primary_field_index("Energy"    )



        !
        ! Get boundary condition input parameters
        !
        coords = worker%coords()
        time   = worker%time()
        rho_input = self%bcproperties%compute("Density",    time, coords)
        p_input   = self%bcproperties%compute("Pressure",   time, coords)
        u_input   = self%bcproperties%compute("X-Velocity", time, coords)
        v_input   = self%bcproperties%compute("Y-Velocity", time, coords)
        w_input   = self%bcproperties%compute("Z-Velocity", time, coords)

        T_input = p_input/(rho_input*287.15_rk)
        c_input = sqrt(1.4_rk*287.15_rk*T_input)



        !
        ! Interpolate interior solution to quadrature nodes
        !
        rho_m  = worker%get_primary_field_face('Density'   ,irho,  'value', 'face interior')
        rhou_m = worker%get_primary_field_face('X-Momentum',irhou, 'value', 'face interior')
        rhov_m = worker%get_primary_field_face('Y-Momentum',irhov, 'value', 'face interior')
        rhow_m = worker%get_primary_field_face('Z-Momentum',irhow, 'value', 'face interior')
        rhoE_m = worker%get_primary_field_face('Energy'    ,irhoE, 'value', 'face interior')

        
        p_m = prop%fluid%compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m)


        T_m = p_m/(rho_m*287.15_rk)
        c_m = sqrt(1.4_rk*287.15_rk*T_m)
!        T_m = prop%fluid%compute_temperature(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m)



        drho_dx_m  = worker%get_primary_field_face('Density'   ,irho,  'ddx', 'face interior')
        drho_dy_m  = worker%get_primary_field_face('Density'   ,irho,  'ddy', 'face interior')
        drho_dz_m  = worker%get_primary_field_face('Density'   ,irho,  'ddz', 'face interior')

        drhou_dx_m = worker%get_primary_field_face('X-Momentum',irhou, 'ddx', 'face interior')
        drhou_dy_m = worker%get_primary_field_face('X-Momentum',irhou, 'ddy', 'face interior')
        drhou_dz_m = worker%get_primary_field_face('X-Momentum',irhou, 'ddz', 'face interior')

        drhov_dx_m = worker%get_primary_field_face('Y-Momentum',irhov, 'ddx', 'face interior')
        drhov_dy_m = worker%get_primary_field_face('Y-Momentum',irhov, 'ddy', 'face interior')
        drhov_dz_m = worker%get_primary_field_face('Y-Momentum',irhov, 'ddz', 'face interior')

        drhow_dx_m = worker%get_primary_field_face('Z-Momentum',irhow, 'ddx', 'face interior')
        drhow_dy_m = worker%get_primary_field_face('Z-Momentum',irhow, 'ddy', 'face interior')
        drhow_dz_m = worker%get_primary_field_face('Z-Momentum',irhow, 'ddz', 'face interior')
        
        drhoE_dx_m = worker%get_primary_field_face('Energy'    ,irhoE, 'ddx', 'face interior')
        drhoE_dy_m = worker%get_primary_field_face('Energy'    ,irhoE, 'ddy', 'face interior')
        drhoE_dz_m = worker%get_primary_field_face('Energy'    ,irhoE, 'ddz', 'face interior')





        ! Initialize arrays
        rho_bc  = rho_m
        rhou_bc = rhou_m
        rhov_bc = rhov_m
        rhow_bc = rhow_m
        rhoE_bc = rhoE_m
        R_inf   = rho_m
        R_extrapolated = rho_m
        u_bc_norm = rho_m
        v_bc_norm = rho_m
        w_bc_norm = rho_m
        u_bc_tang = rho_m
        v_bc_tang = rho_m
        w_bc_tang = rho_m
        entropy_bc = rho_m


        !
        ! Get unit normal vector
        !
        unormx = worker%unit_normal(1)
        unormy = worker%unit_normal(2)
        unormz = worker%unit_normal(3)




        !
        ! Dot momentum with normal vector
        !
        normal_momentum = rhou_m*unormx + rhov_m*unormy + rhow_m*unormz


        !
        ! Determine which nodes are inflow/outflow
        !
        inflow  = ( normal_momentum <= RKTOL )
        outflow = ( normal_momentum >  RKTOL )


        !
        ! Compute internal velocities
        !
        u_m = rhou_m/rho_m
        v_m = rhov_m/rho_m
        w_m = rhow_m/rho_m



        !
        ! Compute Riemann invariants
        !
        R_inf          = (u_input*unormx + v_input*unormy + w_input*unormz) - TWO*c_input/(1.4_rk - ONE)
        R_extrapolated = (u_m*unormx     + v_m*unormy     + w_m*unormz    ) + TWO*c_m/(1.4_rk - ONE)


        !
        ! Compute boundary velocities
        !
        c_bc = ((1.4_rk - ONE)/FOUR)*(R_extrapolated - R_inf)

        u_bc_norm = HALF*(R_extrapolated + R_inf)*unormx
        v_bc_norm = HALF*(R_extrapolated + R_inf)*unormy
        w_bc_norm = HALF*(R_extrapolated + R_inf)*unormz



        !
        ! Compute tangential velocities
        !
        where (inflow)

            u_bc_tang = u_input - (u_input*unormx + v_input*unormy + w_input*unormz)*unormx
            v_bc_tang = v_input - (u_input*unormx + v_input*unormy + w_input*unormz)*unormy
            w_bc_tang = w_input - (u_input*unormx + v_input*unormy + w_input*unormz)*unormz

            entropy_bc = p_input/(rho_input**1.4_rk)

        elsewhere !outflow

            u_bc_tang = u_m - (u_m*unormx + v_m*unormy + w_m*unormz)*unormx
            v_bc_tang = v_m - (u_m*unormx + v_m*unormy + w_m*unormz)*unormy
            w_bc_tang = w_m - (u_m*unormx + v_m*unormy + w_m*unormz)*unormz

            entropy_bc = p_m/(rho_m**1.4_rk)

        end where



        !
        ! Compute boundary state
        !
        rho_bc  = (c_bc*c_bc/(entropy_bc*1.4_rk))**(ONE/(1.4_rk-ONE))
        rhou_bc = (u_bc_norm + u_bc_tang)*rho_bc
        rhov_bc = (v_bc_norm + v_bc_tang)*rho_bc
        rhow_bc = (w_bc_norm + w_bc_tang)*rho_bc

        p_bc   = (rho_bc**1.4_rk)*entropy_bc
        rhoE_bc = (p_bc/(1.4_rk - ONE)) + HALF*(rhou_bc*rhou_bc + rhov_bc*rhov_bc + rhow_bc*rhow_bc)/rho_bc


        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density'   ,irho, rho_bc, 'value')
        call worker%store_bc_state('X-Momentum',irhou,rhou_bc,'value')
        call worker%store_bc_state('Y-Momentum',irhov,rhov_bc,'value')
        call worker%store_bc_state('Z-Momentum',irhow,rhow_bc,'value')
        call worker%store_bc_state('Energy'    ,irhoE,rhoE_bc,'value')


        
        
        
        call worker%store_bc_state('Density'   ,irho,  drho_dx_m,  'ddx')
        call worker%store_bc_state('Density'   ,irho,  drho_dy_m,  'ddy')
        call worker%store_bc_state('Density'   ,irho,  drho_dz_m,  'ddz')
                                                
        call worker%store_bc_state('X-Momentum',irhou, drhou_dx_m, 'ddx')
        call worker%store_bc_state('X-Momentum',irhou, drhou_dy_m, 'ddy')
        call worker%store_bc_state('X-Momentum',irhou, drhou_dz_m, 'ddz')
                                                
        call worker%store_bc_state('Y-Momentum',irhov, drhov_dx_m, 'ddx')
        call worker%store_bc_state('Y-Momentum',irhov, drhov_dy_m, 'ddy')
        call worker%store_bc_state('Y-Momentum',irhov, drhov_dz_m, 'ddz')
                                                
        call worker%store_bc_state('Z-Momentum',irhow, drhow_dx_m, 'ddx')
        call worker%store_bc_state('Z-Momentum',irhow, drhow_dy_m, 'ddy')
        call worker%store_bc_state('Z-Momentum',irhow, drhow_dz_m, 'ddz')
                                                
        call worker%store_bc_state('Energy'    ,irhoE, drhoE_dx_m, 'ddx')
        call worker%store_bc_state('Energy'    ,irhoE, drhoE_dy_m, 'ddy')
        call worker%store_bc_state('Energy'    ,irhoE, drhoE_dz_m, 'ddz')



    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_farfield
