module EULER_boundary_average_advective_flux
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, HALF, ME, NEIGHBOR

    use type_boundary_flux,     only: boundary_flux_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    use EULER_properties,       only: EULER_properties_t
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
    type, extends(boundary_flux_t), public :: EULER_boundary_average_advective_flux_t

    contains

        procedure  :: compute

    end type EULER_boundary_average_advective_flux_t
    !********************************************************************************










contains



    !>  Boundary Flux routine for Euler
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!-------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(EULER_boundary_average_advective_flux_t), intent(in)      :: self
        type(chidg_worker_t),                           intent(inout)   :: worker
        class(properties_t),                            intent(inout)   :: prop

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
            invrho_m,   invrho_p,                   &
            integrand

        real(rk), allocatable, dimension(:) ::      &
            normx, normy, normz


        irho  = prop%get_eqn_index("rho")
        irhou = prop%get_eqn_index("rhou")
        irhov = prop%get_eqn_index("rhov")
        irhow = prop%get_eqn_index("rhow")
        irhoE = prop%get_eqn_index("rhoE")



        !
        ! Interpolate solution to quadrature nodes
        !
        rho_m  = worker%interpolate(irho, 'value', ME)
        rho_p  = worker%interpolate(irho, 'value', NEIGHBOR)

        rhou_m = worker%interpolate(irhou, 'value', ME)
        rhou_p = worker%interpolate(irhou, 'value', NEIGHBOR)

        rhov_m = worker%interpolate(irhov, 'value', ME)
        rhov_p = worker%interpolate(irhov, 'value', NEIGHBOR)

        rhow_m = worker%interpolate(irhow, 'value', ME)
        rhow_p = worker%interpolate(irhow, 'value', NEIGHBOR)

        rhoE_m = worker%interpolate(irhoE, 'value', ME)
        rhoE_p = worker%interpolate(irhoE, 'value', NEIGHBOR)




        invrho_m = ONE/rho_m
        invrho_p = ONE/rho_p



        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)



        !
        ! Compute pressure and total enthalpy
        !
        call prop%fluid%compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,p_m)
        call prop%fluid%compute_pressure(rho_p,rhou_p,rhov_p,rhow_p,rhoE_p,p_p)

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


        ! dot with normal vector
        integrand = HALF*(flux_x*normx + flux_y*normy + flux_z*normz)

        call worker%integrate_boundary(irho, integrand)


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


        ! dot with normal vector
        integrand = HALF*(flux_x*normx + flux_y*normy + flux_z*normz)

        call worker%integrate_boundary(irhou, integrand)


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


        ! dot with normal vector
        integrand = HALF*(flux_x*normx + flux_y*normy + flux_z*normz)

        call worker%integrate_boundary(irhov, integrand)


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


        ! dot with normal vector
        integrand = HALF*(flux_x*normx + flux_y*normy + flux_z*normz)

        call worker%integrate_boundary(irhow, integrand)


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


        ! dot with normal vector
        integrand = HALF*(flux_x*normx + flux_y*normy + flux_z*normz)

        call worker%integrate_boundary(irhoE, integrand)


    end subroutine compute
    !*********************************************************************************************************












end module EULER_boundary_average_advective_flux
