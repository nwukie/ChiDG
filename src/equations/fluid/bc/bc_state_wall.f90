module bc_state_wall
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO
    use type_bc_state,          only: bc_state_t
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
    type, public, extends(bc_state_t) :: wall_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type wall_t
    !*******************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(wall_t),   intent(inout) :: self
        

        !
        ! Set operator name
        !
        call self%set_name("Wall")
        call self%set_family("Wall")


!        !
!        ! Set operator equations
!        !
!        call self%set_equation("Density"   )
!        call self%set_equation("X-Momentum")
!        call self%set_equation("Y-Momentum")
!        call self%set_equation("Z-Momentum")
!        call self%set_equation("Energy"    )


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
        class(wall_t),          intent(inout)   :: self
        type(chidg_worker_t),   intent(inout)   :: worker
        class(properties_t),    intent(inout)   :: prop

        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE

        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::      &
            rho_m,  rhou_m,  rhov_m,  rhow_m,  rhoE_m,  &
            rho_bc, rhou_bc, rhov_bc, rhow_bc, rhoE_bc, &
            drho_dx_m, drhou_dx_m, drhov_dx_m, drhow_dx_m, drhoE_dx_m,  &
            drho_dy_m, drhou_dy_m, drhov_dy_m, drhow_dy_m, drhoE_dy_m,  &
            drho_dz_m, drhou_dz_m, drhov_dz_m, drhow_dz_m, drhoE_dz_m,  &
            u_m, v_m, w_m, VMag2



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





        ! Initialize arrays
        rho_bc  = rho_m
        rhou_bc = rho_m
        rhov_bc = rho_m
        rhow_bc = rho_m
        rhoE_bc = rho_m


        ! Zero momentum
        rhou_bc = ZERO
        rhov_bc = ZERO
        rhow_bc = ZERO


        u_m = rhou_m/rho_m
        v_m = rhov_m/rho_m
        w_m = rhow_m/rho_m

        VMag2 = u_m*u_m + v_m*v_m + w_m*w_m

        !
        ! Energy subtract momentum
        !
        rhoE_bc = rhoE_m - (rho_m*HALF)*VMag2


        !
        ! Store boundary condition state
        !
        call worker%store_bc_state("Density"   ,rho_bc, 'value')
        call worker%store_bc_state("X-Momentum",rhou_bc,'value')
        call worker%store_bc_state("Y-Momentum",rhov_bc,'value')
        call worker%store_bc_state("Z-Momentum",rhow_bc,'value')
        call worker%store_bc_state("Energy"    ,rhoE_bc,'value')


        drho_dx_m = ZERO
        drho_dy_m = ZERO
        drho_dz_m = ZERO
        call worker%store_bc_state("Density"   , drho_dx_m, 'ddx')
        call worker%store_bc_state("Density"   , drho_dy_m, 'ddy')
        call worker%store_bc_state("Density"   , drho_dz_m, 'ddz')
                                                
        call worker%store_bc_state("X-Momentum", drhou_dx_m, 'ddx')
        call worker%store_bc_state("X-Momentum", drhou_dy_m, 'ddy')
        call worker%store_bc_state("X-Momentum", drhou_dz_m, 'ddz')
                                                
        call worker%store_bc_state("Y-Momentum", drhov_dx_m, 'ddx')
        call worker%store_bc_state("Y-Momentum", drhov_dy_m, 'ddy')
        call worker%store_bc_state("Y-Momentum", drhov_dz_m, 'ddz')
                                                
        call worker%store_bc_state("Z-Momentum", drhow_dx_m, 'ddx')
        call worker%store_bc_state("Z-Momentum", drhow_dy_m, 'ddy')
        call worker%store_bc_state("Z-Momentum", drhow_dz_m, 'ddz')

        drhoE_dx_m = ZERO
        drhoE_dy_m = ZERO
        drhoE_dz_m = ZERO
        call worker%store_bc_state("Energy"    , drhoE_dx_m, 'ddx')
        call worker%store_bc_state("Energy"    , drhoE_dy_m, 'ddy')
        call worker%store_bc_state("Energy"    , drhoE_dz_m, 'ddz')



    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_wall
