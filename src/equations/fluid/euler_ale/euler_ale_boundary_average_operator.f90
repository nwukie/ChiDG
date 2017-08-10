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
            invrho_m,   invrho_p,                   &
            integrand

        type(AD_D), allocatable,    dimension(:,:)  ::  &
            flux_ref_m, flux_ref_p

        real(rk), allocatable, dimension(:) ::      &
            normx, normy, normz
        


        !
        ! Interpolate solution to quadrature nodes
        !
        rho_m  = worker%get_primary_field_value_ale_face('Density'   , 'face interior')
        rho_p  = worker%get_primary_field_value_ale_face('Density'   , 'face exterior')

        rhou_m = worker%get_primary_field_value_ale_face('Momentum-1', 'face interior')
        rhou_p = worker%get_primary_field_value_ale_face('Momentum-1', 'face exterior')

        rhov_m = worker%get_primary_field_value_ale_face('Momentum-2', 'face interior')
        rhov_p = worker%get_primary_field_value_ale_face('Momentum-2', 'face exterior')

        rhow_m = worker%get_primary_field_value_ale_face('Momentum-3', 'face interior')
        rhow_p = worker%get_primary_field_value_ale_face('Momentum-3', 'face exterior')

        rhoE_m = worker%get_primary_field_value_ale_face('Energy'    , 'face interior')
        rhoE_p = worker%get_primary_field_value_ale_face('Energy'    , 'face exterior')

        invrho_m = ONE/rho_m
        invrho_p = ONE/rho_p



        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)



        !
        ! Compute pressure and total enthalpy
        !
        p_m = worker%get_model_field_face('Pressure', 'value', 'face interior')
        p_p = worker%get_model_field_face('Pressure', 'value', 'face exterior')


        H_m = (rhoE_m + p_m)*invrho_m
        H_p = (rhoE_p + p_p)*invrho_p



        !================================
        !       MASS FLUX
        !================================
        flux_x_m = rhou_m
        flux_y_m = rhov_m
        flux_z_m = rhow_m
        flux_ref_m = worker%post_process_boundary_advective_flux_ale(flux_x_m,flux_y_m,flux_z_m, advected_quantity=rho_m, interp_source='face interior')

        flux_x_p = rhou_p
        flux_y_p = rhov_p
        flux_z_p = rhow_p
        flux_ref_p = worker%post_process_boundary_advective_flux_ale(flux_x_p,flux_y_p,flux_z_p, advected_quantity=rho_p, interp_source='face exterior')

        flux_x = (flux_ref_m(:,1) + flux_ref_p(:,1))
        flux_y = (flux_ref_m(:,2) + flux_ref_p(:,2))
        flux_z = (flux_ref_m(:,3) + flux_ref_p(:,3))



        ! dot with normal vector
        integrand = HALF*(flux_x*normx + flux_y*normy + flux_z*normz)

        call worker%integrate_boundary('Density',integrand)


        !================================
        !       X-MOMENTUM FLUX
        !================================
        flux_x_m = (rhou_m*rhou_m)*invrho_m + p_m
        flux_y_m = (rhou_m*rhov_m)*invrho_m
        flux_z_m = (rhou_m*rhow_m)*invrho_m
        flux_ref_m = worker%post_process_boundary_advective_flux_ale(flux_x_m,flux_y_m,flux_z_m, advected_quantity=rhou_m, interp_source='face interior')

        flux_x_p = (rhou_p*rhou_p)*invrho_p + p_p
        flux_y_p = (rhou_p*rhov_p)*invrho_p
        flux_z_p = (rhou_p*rhow_p)*invrho_p
        flux_ref_p = worker%post_process_boundary_advective_flux_ale(flux_x_p,flux_y_p,flux_z_p, advected_quantity=rhou_p, interp_source='face exterior')

        flux_x = (flux_ref_m(:,1) + flux_ref_p(:,1))
        flux_y = (flux_ref_m(:,2) + flux_ref_p(:,2))
        flux_z = (flux_ref_m(:,3) + flux_ref_p(:,3))


        ! dot with normal vector
        integrand = HALF*(flux_x*normx + flux_y*normy + flux_z*normz)

        call worker%integrate_boundary('Momentum-1',integrand)


        !================================
        !       Y-MOMENTUM FLUX
        !================================
        flux_x_m = (rhov_m*rhou_m)*invrho_m
        flux_y_m = (rhov_m*rhov_m)*invrho_m + p_m
        flux_z_m = (rhov_m*rhow_m)*invrho_m
        flux_ref_m = worker%post_process_boundary_advective_flux_ale(flux_x_m,flux_y_m,flux_z_m, advected_quantity=rhov_m, interp_source='face interior')

        flux_x_p = (rhov_p*rhou_p)*invrho_p
        flux_y_p = (rhov_p*rhov_p)*invrho_p + p_p
        flux_z_p = (rhov_p*rhow_p)*invrho_p
        flux_ref_p = worker%post_process_boundary_advective_flux_ale(flux_x_p,flux_y_p,flux_z_p, advected_quantity=rhov_p, interp_source='face exterior')

        flux_x = (flux_ref_m(:,1) + flux_ref_p(:,1))
        flux_y = (flux_ref_m(:,2) + flux_ref_p(:,2))
        flux_z = (flux_ref_m(:,3) + flux_ref_p(:,3))



        ! dot with normal vector
        integrand = HALF*(flux_x*normx + flux_y*normy + flux_z*normz)

        call worker%integrate_boundary('Momentum-2',integrand)


        !================================
        !       Z-MOMENTUM FLUX
        !================================
        flux_x_m = (rhow_m*rhou_m)*invrho_m
        flux_y_m = (rhow_m*rhov_m)*invrho_m
        flux_z_m = (rhow_m*rhow_m)*invrho_m + p_m
        flux_ref_m = worker%post_process_boundary_advective_flux_ale(flux_x_m,flux_y_m,flux_z_m, advected_quantity=rhow_m, interp_source='face interior')

        flux_x_p = (rhow_p*rhou_p)*invrho_p
        flux_y_p = (rhow_p*rhov_p)*invrho_p
        flux_z_p = (rhow_p*rhow_p)*invrho_p + p_p
        flux_ref_p = worker%post_process_boundary_advective_flux_ale(flux_x_p,flux_y_p,flux_z_p, advected_quantity=rhow_p, interp_source='face exterior')

        flux_x = (flux_ref_m(:,1) + flux_ref_p(:,1))
        flux_y = (flux_ref_m(:,2) + flux_ref_p(:,2))
        flux_z = (flux_ref_m(:,3) + flux_ref_p(:,3))


        ! dot with normal vector
        integrand = HALF*(flux_x*normx + flux_y*normy + flux_z*normz)

        call worker%integrate_boundary('Momentum-3',integrand)


        !================================
        !          ENERGY FLUX
        !================================
        flux_x_m = rhou_m*H_m
        flux_y_m = rhov_m*H_m
        flux_z_m = rhow_m*H_m
        flux_ref_m = worker%post_process_boundary_advective_flux_ale(flux_x_m,flux_y_m,flux_z_m, advected_quantity=rhoE_m, interp_source='face interior')

        flux_x_p = rhou_p*H_p
        flux_y_p = rhov_p*H_p
        flux_z_p = rhow_p*H_p
        flux_ref_p = worker%post_process_boundary_advective_flux_ale(flux_x_p,flux_y_p,flux_z_p, advected_quantity=rhoE_p, interp_source='face exterior')

        flux_x = (flux_ref_m(:,1) + flux_ref_p(:,1))
        flux_y = (flux_ref_m(:,2) + flux_ref_p(:,2))
        flux_z = (flux_ref_m(:,3) + flux_ref_p(:,3))



        ! dot with normal vector
        integrand = HALF*(flux_x*normx + flux_y*normy + flux_z*normz)

        call worker%integrate_boundary('Energy',integrand)


    end subroutine compute
    !*********************************************************************************************************












end module euler_ale_boundary_average_operator
