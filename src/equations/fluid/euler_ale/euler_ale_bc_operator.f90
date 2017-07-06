module euler_ale_bc_operator
    use mod_constants,      only: HALF,ONE, TWO
    use mod_kinds,          only: ik, rk
    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/12/2016
    !!
    !!
    !------------------------------------------------------------------------------------
    type, public, extends(operator_t)   :: euler_ale_bc_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type euler_ale_bc_operator_t
    !*************************************************************************************






contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_ale_bc_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Euler ALE BC Flux")

        !
        ! Set operator type
        !
        call self%set_operator_type("BC Advective Flux")

        !
        ! Set operator equations
        !
        call self%add_primary_field("Density"   )
        call self%add_primary_field("Momentum-1")
        call self%add_primary_field("Momentum-2")
        call self%add_primary_field("Momentum-3")
        call self%add_primary_field("Energy"    )

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
    subroutine compute(self,worker,prop)
        class(euler_ale_bc_operator_t), intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop


        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::          &
            rho_bc,  rhou_bc, rhov_bc, rhow_bc, rhoE_bc,    &
            u_bc,    v_bc,    w_bc,                         &
            H_bc,    p_bc,                                  &
            flux_x_ref, flux_y_ref, flux_z_ref,     &
            flux_x,  flux_y,  flux_z,  integrand


        real(rk),   allocatable, dimension(:)   ::          &
            normx, normy, normz
            
        real(rk) :: gam_bc

        real(rk), allocatable, dimension(:) ::      &
           u_grid, v_grid, w_grid, det_jacobian_grid, testx


        real(rk), allocatable, dimension(:,:,:) ::      &
            jacobian_grid


        !
        ! Get equation indices
        !
        irho  = prop%get_primary_field_index("Density"   )
        irhou = prop%get_primary_field_index("Momentum-1")
        irhov = prop%get_primary_field_index("Momentum-2")
        irhow = prop%get_primary_field_index("Momentum-3")
        irhoE = prop%get_primary_field_index("Energy"    )



        !
        ! Interpolate boundary condition state to face quadrature nodes
        !
        rho_bc  = worker%get_primary_field_face("Density"   ,'value', 'boundary')
        rhou_bc = worker%get_primary_field_face("Momentum-1",'value', 'boundary')
        rhov_bc = worker%get_primary_field_face("Momentum-2",'value', 'boundary')
        rhow_bc = worker%get_primary_field_face("Momentum-3",'value', 'boundary')
        rhoE_bc = worker%get_primary_field_face("Energy"    ,'value', 'boundary')


        u_grid = worker%get_grid_velocity_face("u_grid")
        v_grid = worker%get_grid_velocity_face("v_grid")
        w_grid = worker%get_grid_velocity_face("w_grid")
        jacobian_grid = worker%get_inv_jacobian_grid_face()
        det_jacobian_grid = worker%get_det_jacobian_grid_face('value')


        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)

        !if (maxval(abs(jacobian_grid(1,:,:))) .ge. 1.002_rk) then
        !print *,'normx'
        !print *,normx(1)
        !print *,'normy'
        !print *,normy(1)
        !print *,'normz'
        !print *,normz(1)

        !print *, 'bc jacobiain'
        !print *, jacobian_grid(1,:,:)
        !end if
        !rho_bc = rho_bc/det_jacobian_grid
        !rhou_bc = rhou_bc/det_jacobian_grid
        !rhov_bc = rhov_bc/det_jacobian_grid
        !rhow_bc = rhow_bc/det_jacobian_grid
        !rhoE_bc = rhoE_bc/det_jacobian_grid



        !
        ! Compute gamma
        !
        !gam_bc = prop%fluid%compute_gamma(rho_bc,rhou_bc,rhov_bc,rhow_bc,rhoE_bc)
        !p_bc   = prop%fluid%compute_pressure(rho_bc,rhou_bc,rhov_bc,rhow_bc,rhoE_bc)
        gam_bc = 1.4_rk
        !p_bc   = worker%get_model_field_face('Pressure','value','boundary')


        p_bc = (1.4_rk-1.0_rk)*(rhoE_bc-HALF*(rhou_bc**TWO+rhov_bc**TWO+rhow_bc**TWO)/rho_bc)



        !
        ! Compute velocity components
        !
        u_bc = rhou_bc/rho_bc
        v_bc = rhov_bc/rho_bc
        w_bc = rhow_bc/rho_bc



        !
        ! Compute boundary condition energy and enthalpy
        !
        H_bc = (rhoE_bc + p_bc)/rho_bc




        !=================================================
        ! Mass flux
        !=================================================
        flux_x = (rho_bc * u_bc) - u_grid*rho_bc 
        flux_y = (rho_bc * v_bc) - v_grid*rho_bc 
        flux_z = (rho_bc * w_bc) - w_grid*rho_bc 

        flux_x_ref = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_x + jacobian_grid(:,1,2)*flux_y + jacobian_grid(:,1,3)*flux_z)
        flux_y_ref = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_x + jacobian_grid(:,2,2)*flux_y + jacobian_grid(:,2,3)*flux_z)
        flux_z_ref = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_x + jacobian_grid(:,3,2)*flux_y + jacobian_grid(:,3,3)*flux_z)


        integrand = flux_x_ref*normx + flux_y_ref*normy + flux_z_ref*normz

        call worker%integrate_boundary('Density',integrand)

        !=================================================
        ! x-momentum flux
        !=================================================
        flux_x = (rho_bc * u_bc * u_bc) - u_grid*rhou_bc  + p_bc
        flux_y = (rho_bc * u_bc * v_bc) - v_grid*rhou_bc 
        flux_z = (rho_bc * u_bc * w_bc) - w_grid*rhou_bc 
        flux_x_ref = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_x + jacobian_grid(:,1,2)*flux_y + jacobian_grid(:,1,3)*flux_z)
        flux_y_ref = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_x + jacobian_grid(:,2,2)*flux_y + jacobian_grid(:,2,3)*flux_z)
        flux_z_ref = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_x + jacobian_grid(:,3,2)*flux_y + jacobian_grid(:,3,3)*flux_z)


        integrand = flux_x_ref*normx + flux_y_ref*normy + flux_z_ref*normz

        call worker%integrate_boundary('Momentum-1',integrand)

        !=================================================
        ! y-momentum flux
        !=================================================
        flux_x = (rho_bc * v_bc * u_bc) - u_grid*rhov_bc 
        flux_y = (rho_bc * v_bc * v_bc) - v_grid*rhov_bc  + p_bc
        flux_z = (rho_bc * v_bc * w_bc) - w_grid*rhov_bc 
        flux_x_ref = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_x + jacobian_grid(:,1,2)*flux_y + jacobian_grid(:,1,3)*flux_z)
        flux_y_ref = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_x + jacobian_grid(:,2,2)*flux_y + jacobian_grid(:,2,3)*flux_z)
        flux_z_ref = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_x + jacobian_grid(:,3,2)*flux_y + jacobian_grid(:,3,3)*flux_z)


        integrand = flux_x_ref*normx + flux_y_ref*normy + flux_z_ref*normz

        call worker%integrate_boundary('Momentum-2',integrand)

        !=================================================
        ! z-momentum flux
        !=================================================
        flux_x = (rho_bc * w_bc * u_bc) - u_grid*rhow_bc 
        flux_y = (rho_bc * w_bc * v_bc) - v_grid*rhow_bc 
        flux_z = (rho_bc * w_bc * w_bc) - w_grid*rhow_bc  + p_bc
        flux_x_ref = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_x + jacobian_grid(:,1,2)*flux_y + jacobian_grid(:,1,3)*flux_z)
        flux_y_ref = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_x + jacobian_grid(:,2,2)*flux_y + jacobian_grid(:,2,3)*flux_z)
        flux_z_ref = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_x + jacobian_grid(:,3,2)*flux_y + jacobian_grid(:,3,3)*flux_z)


        integrand = flux_x_ref*normx + flux_y_ref*normy + flux_z_ref*normz

        call worker%integrate_boundary('Momentum-3',integrand)

        !=================================================
        ! Energy flux
        !=================================================
        flux_x = (rho_bc * u_bc * H_bc) - u_grid*rhoE_bc 
        flux_y = (rho_bc * v_bc * H_bc) - v_grid*rhoE_bc 
        flux_z = (rho_bc * w_bc * H_bc) - w_grid*rhoE_bc 

        flux_x_ref = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_x + jacobian_grid(:,1,2)*flux_y + jacobian_grid(:,1,3)*flux_z)
        flux_y_ref = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_x + jacobian_grid(:,2,2)*flux_y + jacobian_grid(:,2,3)*flux_z)
        flux_z_ref = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_x + jacobian_grid(:,3,2)*flux_y + jacobian_grid(:,3,3)*flux_z)
!        if ((worker%element_info%ielement_g == 1) .and. (worker%iface ==3)) then
!            testx = worker%x('boundary')
!            print *, 'time'
!            print *, worker%t
!            print *, 'node x position'
!            print *, testx(1) 
!            print *, 'det_jacobian_grid'
!            print *, det_jacobian_grid(1)
!            !print *, 'u-grid'
!            !print *, u_grid(1)
!            print *, 'Energy flux sample'
!            print *, flux_x_ref(1)%x_ad_
!        end if



        integrand = flux_x_ref*normx + flux_y_ref*normy + flux_z_ref*normz

        call worker%integrate_boundary('Energy',integrand)

    end subroutine compute
    !**********************************************************************************************























end module euler_ale_bc_operator
