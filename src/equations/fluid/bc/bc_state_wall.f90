module bc_state_wall
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: TWO, HALF, ZERO, ME

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
        rho_m  = worker%get_face_variable(irho,  'value', ME)
        rhou_m = worker%get_face_variable(irhou, 'value', ME)
        rhov_m = worker%get_face_variable(irhov, 'value', ME)
        rhow_m = worker%get_face_variable(irhow, 'value', ME)
        rhoE_m = worker%get_face_variable(irhoE, 'value', ME)


        drho_dx_m  = worker%get_face_variable(irho,  'ddx', ME)
        drho_dy_m  = worker%get_face_variable(irho,  'ddy', ME)
        drho_dz_m  = worker%get_face_variable(irho,  'ddz', ME)

        drhou_dx_m = worker%get_face_variable(irhou, 'ddx', ME)
        drhou_dy_m = worker%get_face_variable(irhou, 'ddy', ME)
        drhou_dz_m = worker%get_face_variable(irhou, 'ddz', ME)

        drhov_dx_m = worker%get_face_variable(irhov, 'ddx', ME)
        drhov_dy_m = worker%get_face_variable(irhov, 'ddy', ME)
        drhov_dz_m = worker%get_face_variable(irhov, 'ddz', ME)

        drhow_dx_m = worker%get_face_variable(irhow, 'ddx', ME)
        drhow_dy_m = worker%get_face_variable(irhow, 'ddy', ME)
        drhow_dz_m = worker%get_face_variable(irhow, 'ddz', ME)
        
        drhoE_dx_m = worker%get_face_variable(irhoE, 'ddx', ME)
        drhoE_dy_m = worker%get_face_variable(irhoE, 'ddy', ME)
        drhoE_dz_m = worker%get_face_variable(irhoE, 'ddz', ME)





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
        call worker%store_bc_state(irho, rho_bc, 'value')
        call worker%store_bc_state(irhou,rhou_bc,'value')
        call worker%store_bc_state(irhov,rhov_bc,'value')
        call worker%store_bc_state(irhow,rhow_bc,'value')
        call worker%store_bc_state(irhoE,rhoE_bc,'value')


        drho_dx_m = ZERO
        drho_dy_m = ZERO
        drho_dz_m = ZERO
        call worker%store_bc_state(irho, drho_dx_m, 'ddx')
        call worker%store_bc_state(irho, drho_dy_m, 'ddy')
        call worker%store_bc_state(irho, drho_dz_m, 'ddz')

        call worker%store_bc_state(irhou, drhou_dx_m, 'ddx')
        call worker%store_bc_state(irhou, drhou_dy_m, 'ddy')
        call worker%store_bc_state(irhou, drhou_dz_m, 'ddz')

        call worker%store_bc_state(irhov, drhov_dx_m, 'ddx')
        call worker%store_bc_state(irhov, drhov_dy_m, 'ddy')
        call worker%store_bc_state(irhov, drhov_dz_m, 'ddz')

        call worker%store_bc_state(irhow, drhow_dx_m, 'ddx')
        call worker%store_bc_state(irhow, drhow_dy_m, 'ddy')
        call worker%store_bc_state(irhow, drhow_dz_m, 'ddz')

        drhoE_dx_m = ZERO
        drhoE_dy_m = ZERO
        drhoE_dz_m = ZERO
        call worker%store_bc_state(irhoE, drhoE_dx_m, 'ddx')
        call worker%store_bc_state(irhoE, drhoE_dy_m, 'ddy')
        call worker%store_bc_state(irhoE, drhoE_dz_m, 'ddz')



    end subroutine compute_bc_state
    !*****************************************************************************************








end module bc_state_wall
