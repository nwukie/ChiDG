module euler_ale_volume_operator
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF

    use type_operator,          only: operator_t
    use type_properties,        only: properties_t
    use type_chidg_worker,      only: chidg_worker_t
    use DNAD_D
    use ieee_arithmetic,        only: ieee_is_nan
    implicit none

    private

    
    !> Volume flux for Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!
    !------------------------------------------------------------------------------
    type, extends(operator_t), public :: euler_ale_volume_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type euler_ale_volume_operator_t
    !******************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(euler_ale_volume_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("Euler ALE Volume Flux")

        ! Set operator type
        call self%set_operator_type("Volume Advective Flux")

        ! Set operator equations
        call self%add_primary_field("Density"   )
        call self%add_primary_field("Momentum-1")
        call self%add_primary_field("Momentum-2")
        call self%add_primary_field("Momentum-3")
        call self%add_primary_field("Energy"    )

    end subroutine init
    !********************************************************************************



    !> Volume flux routine for Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!  
    !!
    !!------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(euler_ale_volume_operator_t), intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop

        ! Equation indices
        integer(ik)    :: irho, irhou, irhov, irhow, irhoE


        type(AD_D), allocatable, dimension(:) ::    &
            rho, rhou, rhov, rhow, rhoE, p, H,      &
            flux_x, flux_y, flux_z, invrho


        type(AD_D), allocatable, dimension(:,:)  :: flux

!        real(rk), allocatable, dimension(:) ::      &
!           u_grid, v_grid, w_grid
!
!
!        real(rk), allocatable, dimension(:,:,:) ::      &
!            jacobian_grid
!
!        u_grid = worker%get_grid_velocity_element("u_grid")
!        v_grid = worker%get_grid_velocity_element("v_grid")
!        w_grid = worker%get_grid_velocity_element("w_grid")
!
!        jacobian_grid = worker%get_inv_jacobian_grid_element()
!        det_jacobian_grid = worker%get_det_jacobian_grid_element('value')


        !
        ! Interpolate solution to quadrature nodes
        !
        rho  = worker%get_primary_field_value_ale_element('Density'   )
        rhou = worker%get_primary_field_value_ale_element('Momentum-1')
        rhov = worker%get_primary_field_value_ale_element('Momentum-2')
        rhow = worker%get_primary_field_value_ale_element('Momentum-3')
        rhoE = worker%get_primary_field_value_ale_element('Energy'    )


        invrho = ONE/rho
    


        !
        ! Compute pressure and total enthalpy
        !
        !p = prop%fluid%compute_pressure(rho,rhou,rhov,rhow,rhoE)
        !p = (1.4_rk-ONE)*(rhoE-HALF*(rhou**TWO+rhov**TWO+rhow**TWO)/rho)
        p = worker%get_model_field_element('Pressure','value')

        H = (rhoE + p)*invrho

        !===========================
        !        MASS FLUX
        !===========================
        flux_x = rhou
        flux_y = rhov
        flux_z = rhow
        flux = worker%post_process_volume_advective_flux_ale(flux_x,flux_y,flux_z, advected_quantity=rho)

!        flux_x = flux_x - rho*u_grid
!        flux_y = flux_y - rho*v_grid
!        flux_z = flux_z - rho*w_grid
!
!        flux_x_ref = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_x + jacobian_grid(:,1,2)*flux_y + jacobian_grid(:,1,3)*flux_z)
!        flux_y_ref = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_x + jacobian_grid(:,2,2)*flux_y + jacobian_grid(:,2,3)*flux_z)
!        flux_z_ref = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_x + jacobian_grid(:,3,2)*flux_y + jacobian_grid(:,3,3)*flux_z)


        !call worker%integrate_volume('Density',flux_x_ref,flux_y_ref,flux_z_ref)
        call worker%integrate_volume('Density',flux(:,1),flux(:,2),flux(:,3))


        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux_x = (rhou*rhou)*invrho  +  p
        flux_y = (rhou*rhov)*invrho
        flux_z = (rhou*rhow)*invrho
        flux = worker%post_process_volume_advective_flux_ale(flux_x,flux_y,flux_z, advected_quantity=rhou)
        
!        flux_x = flux_x - rhou*u_grid
!        flux_y = flux_y - rhou*v_grid
!        flux_z = flux_z - rhou*w_grid
!
!        flux_x_ref = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_x + jacobian_grid(:,1,2)*flux_y + jacobian_grid(:,1,3)*flux_z)
!        flux_y_ref = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_x + jacobian_grid(:,2,2)*flux_y + jacobian_grid(:,2,3)*flux_z)
!        flux_z_ref = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_x + jacobian_grid(:,3,2)*flux_y + jacobian_grid(:,3,3)*flux_z)


        !call worker%integrate_volume('Momentum-1',flux_x_ref,flux_y_ref,flux_z_ref)
        call worker%integrate_volume('Momentum-1',flux(:,1),flux(:,2),flux(:,3))


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux_x = (rhov*rhou)*invrho
        flux_y = (rhov*rhov)*invrho  +  p
        flux_z = (rhov*rhow)*invrho
        flux = worker%post_process_volume_advective_flux_ale(flux_x,flux_y,flux_z, advected_quantity=rhov)
        
!        flux_x = flux_x - rhov*u_grid
!        flux_y = flux_y - rhov*v_grid
!        flux_z = flux_z - rhov*w_grid
!
!        flux_x_ref = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_x + jacobian_grid(:,1,2)*flux_y + jacobian_grid(:,1,3)*flux_z)
!        flux_y_ref = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_x + jacobian_grid(:,2,2)*flux_y + jacobian_grid(:,2,3)*flux_z)
!        flux_z_ref = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_x + jacobian_grid(:,3,2)*flux_y + jacobian_grid(:,3,3)*flux_z)


        !call worker%integrate_volume('Momentum-2',flux_x_ref,flux_y_ref,flux_z_ref)
        call worker%integrate_volume('Momentum-2',flux(:,1),flux(:,2),flux(:,3))

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        flux_x = (rhow*rhou)*invrho
        flux_y = (rhow*rhov)*invrho
        flux_z = (rhow*rhow)*invrho  +  p
        flux = worker%post_process_volume_advective_flux_ale(flux_x,flux_y,flux_z, advected_quantity=rhow)

!        flux_x = flux_x - rhow*u_grid
!        flux_y = flux_y - rhow*v_grid
!        flux_z = flux_z - rhow*w_grid
!
!        flux_x_ref = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_x + jacobian_grid(:,1,2)*flux_y + jacobian_grid(:,1,3)*flux_z)
!        flux_y_ref = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_x + jacobian_grid(:,2,2)*flux_y + jacobian_grid(:,2,3)*flux_z)
!        flux_z_ref = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_x + jacobian_grid(:,3,2)*flux_y + jacobian_grid(:,3,3)*flux_z)


        !call worker%integrate_volume('Momentum-3',flux_x_ref,flux_y_ref,flux_z_ref)
        call worker%integrate_volume('Momentum-3',flux(:,1),flux(:,2),flux(:,3))

        !============================
        !       ENERGY FLUX
        !============================
        flux_x = rhou*H
        flux_y = rhov*H
        flux_z = rhow*H
        flux = worker%post_process_volume_advective_flux_ale(flux_x,flux_y,flux_z, advected_quantity=rhoE)

!        flux_x = flux_x - rhoE*u_grid
!        flux_y = flux_y - rhoE*v_grid
!        flux_z = flux_z - rhoE*w_grid
!
!        flux_x_ref = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_x + jacobian_grid(:,1,2)*flux_y + jacobian_grid(:,1,3)*flux_z)
!        flux_y_ref = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_x + jacobian_grid(:,2,2)*flux_y + jacobian_grid(:,2,3)*flux_z)
!        flux_z_ref = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_x + jacobian_grid(:,3,2)*flux_y + jacobian_grid(:,3,3)*flux_z)


        !call worker%integrate_volume('Energy',flux_x_ref,flux_y_ref,flux_z_ref)
        call worker%integrate_volume('Energy',flux(:,1),flux(:,2),flux(:,3))

    end subroutine compute
    !*********************************************************************************************************






end module euler_ale_volume_operator
