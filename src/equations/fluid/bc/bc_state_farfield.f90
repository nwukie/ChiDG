module bc_state_farfield
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO, ONE, RKTOL, FOUR

    use type_bc_state,          only: bc_state_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use type_point,             only: point_t
    use mpi_f08,                only: mpi_comm
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



        !
        ! Add boundary condition parameters
        !
        call self%bcproperties%add('Density',   'Required')
        call self%bcproperties%add('Pressure',  'Required')
        call self%bcproperties%add('Velocity-1','Required')
        call self%bcproperties%add('Velocity-2','Required')
        call self%bcproperties%add('Velocity-3','Required')



        !
        ! Set default parameter values
        !
        call self%set_fcn_option('Density',    'val', 1.2_rk)
        call self%set_fcn_option('Pressure',   'val', 100000._rk)
        call self%set_fcn_option('Velocity-1', 'val', 0._rk)
        call self%set_fcn_option('Velocity-2', 'val', 0._rk)
        call self%set_fcn_option('Velocity-3', 'val', 0._rk)



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
    subroutine compute_bc_state(self,worker,prop,bc_COMM)
        class(farfield_t),      intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        class(properties_t),    intent(inout)   :: prop
        type(mpi_comm),         intent(in)      :: bc_COMM

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::                      &
            density_m,  mom1_m,  mom2_m,  mom3_m,  energy_m,            &
            density_bc, mom1_bc, mom2_bc, mom3_bc, energy_bc, p_bc,     &
            drho_dx_m, drhou_dx_m, drhov_dx_m, drhow_dx_m, drhoE_dx_m,  &
            drho_dy_m, drhou_dy_m, drhov_dy_m, drhow_dy_m, drhoE_dy_m,  &
            drho_dz_m, drhou_dz_m, drhov_dz_m, drhow_dz_m, drhoE_dz_m,  &
            normal_momentum, normal_velocity, R_inf, R_extrapolated,    &
            u_bc_norm, v_bc_norm, w_bc_norm, u_bc_tang, v_bc_tang, w_bc_tang, entropy_bc, &
            c_bc, c_m, p_m, T_m, u_m, v_m, w_m

        real(rk), allocatable, dimension(:) ::              &
            unorm_1, unorm_2, unorm_3, r,                   &
            rho_input, p_input, u_input, v_input, w_input, T_input, c_input

        real(rk), allocatable, dimension(:) ::      &
           u_grid, v_grid, w_grid, det_jacobian_grid

        logical, allocatable, dimension(:)  :: inflow, outflow


        det_jacobian_grid = worker%get_det_jacobian_grid_face()
        !
        ! Get boundary condition input parameters
        !
        rho_input = self%bcproperties%compute('Density',    worker%time(), worker%coords())
        p_input   = self%bcproperties%compute('Pressure',   worker%time(), worker%coords())
        u_input   = self%bcproperties%compute('Velocity-1', worker%time(), worker%coords())
        v_input   = self%bcproperties%compute('Velocity-2', worker%time(), worker%coords())
        w_input   = self%bcproperties%compute('Velocity-3', worker%time(), worker%coords())

        T_input = p_input/(rho_input*287.15_rk)
        c_input = sqrt(1.4_rk*287.15_rk*T_input)



        !
        ! Interpolate interior solution to quadrature nodes
        !
        density_m = worker%get_primary_field_face('Density'   , 'value', 'face interior')
        mom1_m    = worker%get_primary_field_face('Momentum-1', 'value', 'face interior')
        mom2_m    = worker%get_primary_field_face('Momentum-2', 'value', 'face interior')
        mom3_m    = worker%get_primary_field_face('Momentum-3', 'value', 'face interior')
        energy_m  = worker%get_primary_field_face('Energy'    , 'value', 'face interior')

        density_m = density_m/det_jacobian_grid
        mom1_m = mom1_m/det_jacobian_grid
        mom2_m = mom2_m/det_jacobian_grid
        mom3_m = mom3_m/det_jacobian_grid
        energy_m = energy_m/det_jacobian_grid



        !
        ! Account for cylindrical. Get tangential momentum from angular momentum.
        !
        r = worker%coordinate('1','boundary')
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_m = mom2_m / r
        end if



        !
        ! Get Pressure, Temperature from interior
        !
        p_m = worker%get_model_field_face('Pressure',    'value', 'face interior')
        T_m = worker%get_model_field_face('Temperature', 'value', 'face interior')


        c_m = sqrt(1.4_rk*287.15_rk*T_m)



        drho_dx_m  = worker%get_primary_field_face('Density'   , 'grad1', 'face interior')
        drho_dy_m  = worker%get_primary_field_face('Density'   , 'grad2', 'face interior')
        drho_dz_m  = worker%get_primary_field_face('Density'   , 'grad3', 'face interior')

        drhou_dx_m = worker%get_primary_field_face('Momentum-1', 'grad1', 'face interior')
        drhou_dy_m = worker%get_primary_field_face('Momentum-1', 'grad2', 'face interior')
        drhou_dz_m = worker%get_primary_field_face('Momentum-1', 'grad3', 'face interior')

        drhov_dx_m = worker%get_primary_field_face('Momentum-2', 'grad1', 'face interior')
        drhov_dy_m = worker%get_primary_field_face('Momentum-2', 'grad2', 'face interior')
        drhov_dz_m = worker%get_primary_field_face('Momentum-2', 'grad3', 'face interior')

        drhow_dx_m = worker%get_primary_field_face('Momentum-3', 'grad1', 'face interior')
        drhow_dy_m = worker%get_primary_field_face('Momentum-3', 'grad2', 'face interior')
        drhow_dz_m = worker%get_primary_field_face('Momentum-3', 'grad3', 'face interior')
        
        drhoE_dx_m = worker%get_primary_field_face('Energy'    , 'grad1', 'face interior')
        drhoE_dy_m = worker%get_primary_field_face('Energy'    , 'grad2', 'face interior')
        drhoE_dz_m = worker%get_primary_field_face('Energy'    , 'grad3', 'face interior')





        ! Initialize arrays
        density_bc = density_m
        mom1_bc    = mom1_m
        mom2_bc    = mom2_m
        mom3_bc    = mom3_m
        energy_bc  = energy_m
        R_inf      = density_m
        R_extrapolated = density_m
        u_bc_norm = density_m
        v_bc_norm = density_m
        w_bc_norm = density_m
        u_bc_tang = density_m
        v_bc_tang = density_m
        w_bc_tang = density_m
        entropy_bc = density_m


        !
        ! Get unit normal vector
        !
        unorm_1 = worker%unit_normal(1)
        unorm_2 = worker%unit_normal(2)
        unorm_3 = worker%unit_normal(3)




        !
        ! Dot momentum with normal vector
        !
        normal_momentum = mom1_m*unorm_1 + mom2_m*unorm_2 + mom3_m*unorm_3


        !
        ! Determine which nodes are inflow/outflow
        !
        inflow  = ( normal_momentum <= RKTOL )
        outflow = ( normal_momentum >  RKTOL )


        !
        ! Compute internal velocities
        !
        u_m = mom1_m/density_m
        v_m = mom2_m/density_m
        w_m = mom3_m/density_m



        !
        ! Compute Riemann invariants
        !
        R_inf          = (u_input*unorm_1 + v_input*unorm_2 + w_input*unorm_3) - TWO*c_input/(1.4_rk - ONE)
        R_extrapolated = (u_m*unorm_1     + v_m*unorm_2     + w_m*unorm_3    ) + TWO*c_m/(1.4_rk - ONE)


        !
        ! Compute boundary velocities
        !
        c_bc = ((1.4_rk - ONE)/FOUR)*(R_extrapolated - R_inf)

        u_bc_norm = HALF*(R_extrapolated + R_inf)*unorm_1
        v_bc_norm = HALF*(R_extrapolated + R_inf)*unorm_2
        w_bc_norm = HALF*(R_extrapolated + R_inf)*unorm_3



        !
        ! Compute tangential velocities
        !
        where (inflow)

            u_bc_tang = u_input - (u_input*unorm_1 + v_input*unorm_2 + w_input*unorm_3)*unorm_1
            v_bc_tang = v_input - (u_input*unorm_1 + v_input*unorm_2 + w_input*unorm_3)*unorm_2
            w_bc_tang = w_input - (u_input*unorm_1 + v_input*unorm_2 + w_input*unorm_3)*unorm_3

            entropy_bc = p_input/(rho_input**1.4_rk)

        elsewhere !outflow

            u_bc_tang = u_m - (u_m*unorm_1 + v_m*unorm_2 + w_m*unorm_3)*unorm_1
            v_bc_tang = v_m - (u_m*unorm_1 + v_m*unorm_2 + w_m*unorm_3)*unorm_2
            w_bc_tang = w_m - (u_m*unorm_1 + v_m*unorm_2 + w_m*unorm_3)*unorm_3

            entropy_bc = p_m/(density_m**1.4_rk)

        end where



        !
        ! Compute boundary state
        !
        density_bc  = (c_bc*c_bc/(entropy_bc*1.4_rk))**(ONE/(1.4_rk-ONE))
        mom1_bc = (u_bc_norm + u_bc_tang)*density_bc
        mom2_bc = (v_bc_norm + v_bc_tang)*density_bc
        mom3_bc = (w_bc_norm + w_bc_tang)*density_bc

        p_bc   = (density_bc**1.4_rk)*entropy_bc
        energy_bc = (p_bc/(1.4_rk - ONE)) + HALF*(mom1_bc*mom1_bc + mom2_bc*mom2_bc + mom3_bc*mom3_bc)/density_bc




        !
        ! Account for cylindrical. Convert tangential momentum back to angular momentum.
        !
        if (worker%coordinate_system() == 'Cylindrical') then
            mom2_bc = mom2_bc * r
        end if


        !
        ! Store boundary condition state
        !
        call worker%store_bc_state('Density'   , density_bc, 'value')
        call worker%store_bc_state('Momentum-1', mom1_bc,    'value')
        call worker%store_bc_state('Momentum-2', mom2_bc,    'value')
        call worker%store_bc_state('Momentum-3', mom3_bc,    'value')
        call worker%store_bc_state('Energy'    , energy_bc,  'value')


        
        
        
        call worker%store_bc_state('Density'   , drho_dx_m,  'grad1')
        call worker%store_bc_state('Density'   , drho_dy_m,  'grad2')
        call worker%store_bc_state('Density'   , drho_dz_m,  'grad3')
                                                
        call worker%store_bc_state('Momentum-1', drhou_dx_m, 'grad1')
        call worker%store_bc_state('Momentum-1', drhou_dy_m, 'grad2')
        call worker%store_bc_state('Momentum-1', drhou_dz_m, 'grad3')
                                                
        call worker%store_bc_state('Momentum-2', drhov_dx_m, 'grad1')
        call worker%store_bc_state('Momentum-2', drhov_dy_m, 'grad2')
        call worker%store_bc_state('Momentum-2', drhov_dz_m, 'grad3')
                                                
        call worker%store_bc_state('Momentum-3', drhow_dx_m, 'grad1')
        call worker%store_bc_state('Momentum-3', drhow_dy_m, 'grad2')
        call worker%store_bc_state('Momentum-3', drhow_dz_m, 'grad3')
                                                
        call worker%store_bc_state('Energy'    , drhoE_dx_m, 'grad1')
        call worker%store_bc_state('Energy'    , drhoE_dy_m, 'grad2')
        call worker%store_bc_state('Energy'    , drhoE_dz_m, 'grad3')



    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_farfield
