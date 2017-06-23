module euler_ale_boundary_average_operator
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, HALF
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none

    private



    !> Implementation of the Euler boundary average flux
    !!
    !!  - At a boundary interface, the solution states Q- and Q+ exists on opposite 
    !!    sides of the boundary. The average flux is computed as Favg = 1/2(F(Q-) + F(Q+))
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: euler_ale_boundary_average_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type euler_ale_boundary_average_operator_t
    !********************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_ale_boundary_average_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Euler ALE Boundary Average Flux")

        !
        ! Set operator type
        !
        call self%set_operator_type("Boundary Advective Flux")

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



    !>  Boundary Flux routine for Euler
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!-------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(euler_ale_boundary_average_operator_t),   intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop

        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE


        ! Storage at quadrature nodes
        type(AD_D), allocatable,    dimension(:) :: &
            rho_m,      rho_p,                      &
            rhou_m,     rhou_p,                     &
            rhov_m,     rhov_p,                     &
            rhow_m,     rhow_p,                     &
            rhoe_m,     rhoe_p,                     &
            rhor_p,     rhot_p,                     &
            p_m,        p_p,                        &
            H_m,        H_p,                        &
            flux_x_m,   flux_y_m,   flux_z_m,       &
            flux_x_p,   flux_y_p,   flux_z_p,       &
            flux_x,     flux_y,     flux_z,         &
            flux_x_ref, flux_y_ref, flux_z_ref,     &
            invrho_m,   invrho_p,                   &
            integrand

        real(rk), allocatable, dimension(:) ::      &
            normx, normy, normz
        
        real(rk), allocatable, dimension(:) ::      &
           u_grid, v_grid, w_grid, det_jacobian_grid


        real(rk), allocatable, dimension(:,:,:) ::      &
            jacobian_grid

        u_grid = worker%get_grid_velocity_face("u_grid")
        v_grid = worker%get_grid_velocity_face("v_grid")
        w_grid = worker%get_grid_velocity_face("w_grid")

        jacobian_grid = worker%get_inv_jacobian_grid_face()
        det_jacobian_grid = worker%get_det_jacobian_grid_face()
!        print *, 'jacobian_grid'
!        print *, jacobian_grid(2,:,:)
!        print *, 'det_jacobian_grid'
!        print *, det_jacobian_grid


        ! Or do we want a single procedure?
        ! worker%get_ale_field_face("field_name")

        irho  = prop%get_primary_field_index("Density"   )
        irhou = prop%get_primary_field_index("Momentum-1")
        irhov = prop%get_primary_field_index("Momentum-2")
        irhow = prop%get_primary_field_index("Momentum-3")
        irhoE = prop%get_primary_field_index("Energy"    )



        !
        ! Interpolate solution to quadrature nodes
        !
        rho_m  = worker%get_primary_field_face('Density'   , 'value', 'face interior')
        rho_p  = worker%get_primary_field_face('Density'   , 'value', 'face exterior')

        rhou_m = worker%get_primary_field_face('Momentum-1', 'value', 'face interior')
        rhou_p = worker%get_primary_field_face('Momentum-1', 'value', 'face exterior')

        rhov_m = worker%get_primary_field_face('Momentum-2', 'value', 'face interior')
        rhov_p = worker%get_primary_field_face('Momentum-2', 'value', 'face exterior')

        rhow_m = worker%get_primary_field_face('Momentum-3', 'value', 'face interior')
        rhow_p = worker%get_primary_field_face('Momentum-3', 'value', 'face exterior')

        rhoE_m = worker%get_primary_field_face('Energy'    , 'value', 'face interior')
        rhoE_p = worker%get_primary_field_face('Energy'    , 'value', 'face exterior')

        rho_m = rho_m/det_jacobian_grid
        rhou_m = rhou_m/det_jacobian_grid
        rhov_m = rhov_m/det_jacobian_grid
        rhow_m = rhow_m/det_jacobian_grid
        rhoE_m = rhoE_m/det_jacobian_grid

        rho_p = rho_p/det_jacobian_grid
        rhou_p = rhou_p/det_jacobian_grid
        rhov_p = rhov_p/det_jacobian_grid
        rhow_p = rhow_p/det_jacobian_grid
        rhoE_p = rhoE_p/det_jacobian_grid

        invrho_m = ONE/rho_m
        invrho_p = ONE/rho_p



        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)



        !
        ! Compute pressure and total enthalpy
        !
        !p_m = prop%fluid%compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m)
        !p_p = prop%fluid%compute_pressure(rho_p,rhou_p,rhov_p,rhow_p,rhoE_p)
        !p_m = worker%get_model_field_face('Pressure', 'value', 'face interior')
        !p_p = worker%get_model_field_face('Pressure', 'value', 'face exterior')

        p_m = (1.4_rk-ONE)*(rhoE_m - HALF*(rhou_m**TWO+rhov_m**TWO+rhow_m**TWO)/rho_m)
        p_p = (1.4_rk-ONE)*(rhoE_p - HALF*(rhou_p**TWO+rhov_p**TWO+rhow_p**TWO)/rho_p)

        H_m = (rhoE_m + p_m)*invrho_m
        H_p = (rhoE_p + p_p)*invrho_p



        !================================
        !       MASS FLUX
        !================================
        flux_x_m = rhou_m
        flux_y_m = rhov_m
        flux_z_m = rhow_m

        flux_x_p = rhou_p
        flux_y_p = rhov_p
        flux_z_p = rhow_p

        flux_x = (flux_x_m + flux_x_p)
        flux_y = (flux_y_m + flux_y_p)
        flux_z = (flux_z_m + flux_z_p)
        
        flux_x = flux_x - (rho_m+rho_p)*u_grid
        flux_y = flux_y - (rho_m+rho_p)*v_grid
        flux_z = flux_z - (rho_m+rho_p)*w_grid

        flux_x_ref = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_x + jacobian_grid(:,1,2)*flux_y + jacobian_grid(:,1,3)*flux_z)
        flux_y_ref = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_x + jacobian_grid(:,2,2)*flux_y + jacobian_grid(:,2,3)*flux_z)
        flux_z_ref = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_x + jacobian_grid(:,3,2)*flux_y + jacobian_grid(:,3,3)*flux_z)



        ! dot with normal vector
        integrand = HALF*(flux_x_ref*normx + flux_y_ref*normy + flux_z_ref*normz)

        call worker%integrate_boundary('Density',integrand)


        !================================
        !       X-MOMENTUM FLUX
        !================================
        flux_x_m = (rhou_m*rhou_m)*invrho_m + p_m
        flux_y_m = (rhou_m*rhov_m)*invrho_m
        flux_z_m = (rhou_m*rhow_m)*invrho_m

        flux_x_p = (rhou_p*rhou_p)*invrho_p + p_p
        flux_y_p = (rhou_p*rhov_p)*invrho_p
        flux_z_p = (rhou_p*rhow_p)*invrho_p

        flux_x = (flux_x_m + flux_x_p)
        flux_y = (flux_y_m + flux_y_p)
        flux_z = (flux_z_m + flux_z_p)
        
        flux_x = flux_x - (rhou_m+rhou_p)*u_grid
        flux_y = flux_y - (rhou_m+rhou_p)*v_grid
        flux_z = flux_z - (rhou_m+rhou_p)*w_grid

        flux_x_ref = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_x + jacobian_grid(:,1,2)*flux_y + jacobian_grid(:,1,3)*flux_z)
        flux_y_ref = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_x + jacobian_grid(:,2,2)*flux_y + jacobian_grid(:,2,3)*flux_z)
        flux_z_ref = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_x + jacobian_grid(:,3,2)*flux_y + jacobian_grid(:,3,3)*flux_z)



        ! dot with normal vector
        integrand = HALF*(flux_x_ref*normx + flux_y_ref*normy + flux_z_ref*normz)

        call worker%integrate_boundary('Momentum-1',integrand)


        !================================
        !       Y-MOMENTUM FLUX
        !================================
        flux_x_m = (rhov_m*rhou_m)*invrho_m
        flux_y_m = (rhov_m*rhov_m)*invrho_m + p_m
        flux_z_m = (rhov_m*rhow_m)*invrho_m

        flux_x_p = (rhov_p*rhou_p)*invrho_p
        flux_y_p = (rhov_p*rhov_p)*invrho_p + p_p
        flux_z_p = (rhov_p*rhow_p)*invrho_p

        flux_x = (flux_x_m + flux_x_p)
        flux_y = (flux_y_m + flux_y_p)
        flux_z = (flux_z_m + flux_z_p)

        flux_x = flux_x - (rhov_m+rhov_p)*u_grid
        flux_y = flux_y - (rhov_m+rhov_p)*v_grid
        flux_z = flux_z - (rhov_m+rhov_p)*w_grid

        flux_x_ref = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_x + jacobian_grid(:,1,2)*flux_y + jacobian_grid(:,1,3)*flux_z)
        flux_y_ref = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_x + jacobian_grid(:,2,2)*flux_y + jacobian_grid(:,2,3)*flux_z)
        flux_z_ref = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_x + jacobian_grid(:,3,2)*flux_y + jacobian_grid(:,3,3)*flux_z)


        ! dot with normal vector
        integrand = HALF*(flux_x_ref*normx + flux_y_ref*normy + flux_z_ref*normz)

        call worker%integrate_boundary('Momentum-2',integrand)


        !================================
        !       Z-MOMENTUM FLUX
        !================================
        flux_x_m = (rhow_m*rhou_m)*invrho_m
        flux_y_m = (rhow_m*rhov_m)*invrho_m
        flux_z_m = (rhow_m*rhow_m)*invrho_m + p_m

        flux_x_p = (rhow_p*rhou_p)*invrho_p
        flux_y_p = (rhow_p*rhov_p)*invrho_p
        flux_z_p = (rhow_p*rhow_p)*invrho_p + p_p

        flux_x = (flux_x_m + flux_x_p)
        flux_y = (flux_y_m + flux_y_p)
        flux_z = (flux_z_m + flux_z_p)

        flux_x = flux_x - (rhow_m+rhow_p)*u_grid
        flux_y = flux_y - (rhow_m+rhow_p)*v_grid
        flux_z = flux_z - (rhow_m+rhow_p)*w_grid

        flux_x_ref = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_x + jacobian_grid(:,1,2)*flux_y + jacobian_grid(:,1,3)*flux_z)
        flux_y_ref = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_x + jacobian_grid(:,2,2)*flux_y + jacobian_grid(:,2,3)*flux_z)
        flux_z_ref = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_x + jacobian_grid(:,3,2)*flux_y + jacobian_grid(:,3,3)*flux_z)


        ! dot with normal vector
        integrand = HALF*(flux_x_ref*normx + flux_y_ref*normy + flux_z_ref*normz)

        call worker%integrate_boundary('Momentum-3',integrand)


        !================================
        !          ENERGY FLUX
        !================================
        flux_x_m = rhou_m*H_m
        flux_y_m = rhov_m*H_m
        flux_z_m = rhow_m*H_m

        flux_x_p = rhou_p*H_p
        flux_y_p = rhov_p*H_p
        flux_z_p = rhow_p*H_p

        flux_x = (flux_x_m + flux_x_p)
        flux_y = (flux_y_m + flux_y_p)
        flux_z = (flux_z_m + flux_z_p)

        flux_x = flux_x - (rhoE_m+rhoE_p)*u_grid
        flux_y = flux_y - (rhoE_m+rhoE_p)*v_grid
        flux_z = flux_z - (rhoE_m+rhoE_p)*w_grid

        flux_x_ref = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_x + jacobian_grid(:,1,2)*flux_y + jacobian_grid(:,1,3)*flux_z)
        flux_y_ref = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_x + jacobian_grid(:,2,2)*flux_y + jacobian_grid(:,2,3)*flux_z)
        flux_z_ref = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_x + jacobian_grid(:,3,2)*flux_y + jacobian_grid(:,3,3)*flux_z)


        ! dot with normal vector
        integrand = HALF*(flux_x_ref*normx + flux_y_ref*normy + flux_z_ref*normz)

        call worker%integrate_boundary('Energy',integrand)


    end subroutine compute
    !*********************************************************************************************************












end module euler_ale_boundary_average_operator
