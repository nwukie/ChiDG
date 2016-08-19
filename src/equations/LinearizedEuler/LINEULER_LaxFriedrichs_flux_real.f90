module LINEULER_LaxFriedrichs_flux_real
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF,FOUR,ZERO, ME, NEIGHBOR

    use type_boundary_flux,     only: boundary_flux_t
    use type_mesh,              only: mesh_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    use type_face_info,         only: face_info_t
    use type_function_info,     only: function_info_t

    use mod_interpolate,        only: interpolate
    use mod_integrate,          only: integrate_boundary_scalar_flux
    use DNAD_D

    use LINEULER_properties,    only: LINEULER_properties_t
    use mod_linearized_euler
    implicit none

    private











    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!
    !-------------------------------------------------------------------------------------
    type, extends(boundary_flux_t), public :: LINEULER_LaxFriedrichs_flux_real_t

    contains

        procedure  :: compute

    end type LINEULER_LaxFriedrichs_flux_real_t
    !**************************************************************************************










contains






    !>  Real component of numerical flux dissipation
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,face_info,function_info)
        class(LINEULER_LaxFriedrichs_flux_real_t),  intent(in)      :: self
        type(mesh_t),                               intent(in)      :: mesh(:)
        type(solverdata_t),                         intent(inout)   :: sdata
        class(properties_t),                        intent(inout)   :: prop
        type(face_info_t),                          intent(in)      :: face_info
        type(function_info_t),                      intent(in)      :: function_info

        ! Equation indices
        integer(ik)     :: irho
        integer(ik)     :: irhou
        integer(ik)     :: irhov
        integer(ik)     :: irhow
        integer(ik)     :: irhoe

        integer(ik)     :: idom, ielem,  iface, igq


        real(rk)    :: a_c


        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(face_info%idomain_l)%faces(face_info%ielement_l,face_info%iface)%gq%face%nnodes)    :: &
                        rho_m,      rho_p,                                        &
                        rhou_m,     rhou_p,                                       &
                        rhov_m,     rhov_p,                                       &
                        rhow_m,     rhow_p,                                       &
                        rhoe_m,     rhoe_p,                                       &
                        integrand,  upwind,                                       &
                        wave


        real(rk), dimension(mesh(face_info%idomain_l)%faces(face_info%ielement_l,face_info%iface)%gq%face%nnodes)    :: &
                        un, wave_c


        irho  = prop%get_eqn_index("rho_r")
        irhou = prop%get_eqn_index("rhou_r")
        irhov = prop%get_eqn_index("rhov_r")
        irhow = prop%get_eqn_index("rhow_r")
        irhoE = prop%get_eqn_index("rhoE_r")


        idom  = face_info%idomain_l
        ielem = face_info%ielement_l
        iface = face_info%iface


        associate (norms => mesh(idom)%faces(ielem,iface)%norm, unorms=> mesh(idom)%faces(ielem,iface)%unorm, faces => mesh(idom)%faces, q => sdata%q)

            !
            ! Interpolate solution to quadrature nodes
            !
            rho_m  = interpolate(mesh,sdata,face_info,function_info, irho,  'value',  ME)
            rho_p  = interpolate(mesh,sdata,face_info,function_info, irho,  'value',  NEIGHBOR)

            rhou_m = interpolate(mesh,sdata,face_info,function_info, irhou, 'value', ME)
            rhou_p = interpolate(mesh,sdata,face_info,function_info, irhou, 'value', NEIGHBOR)

            rhov_m = interpolate(mesh,sdata,face_info,function_info, irhov, 'value', ME)
            rhov_p = interpolate(mesh,sdata,face_info,function_info, irhov, 'value', NEIGHBOR)

            rhow_m = interpolate(mesh,sdata,face_info,function_info, irhow, 'value', ME)
            rhow_p = interpolate(mesh,sdata,face_info,function_info, irhow, 'value', NEIGHBOR)

            rhoE_m = interpolate(mesh,sdata,face_info,function_info, irhoE, 'value', ME)
            rhoE_p = interpolate(mesh,sdata,face_info,function_info, irhoE, 'value', NEIGHBOR)

            wave = rho_m
            do igq = 1,size(wave)
                wave(igq)%x_ad_ = ZERO
                wave(igq)%xp_ad_ = ZERO
            end do

            !--------------------------------------
            !  Compute wave speeds
            !--------------------------------------

            !
            ! Compute speed of sound
            !
            a_c = sqrt(abs(gam * pbar / rho_c))


            !
            ! Compute normal velocities: dot-product vector projection along unit-normal direction
            !
            un = unorms(:,1)*ubar + unorms(:,2)*vbar + unorms(:,3)*wbar

            !
            ! Compute wave speeds
            !
            wave_c = abs(un) + a_c

            wave = wave_c


            !================================
            !       MASS FLUX
            !================================
            upwind = -wave*(rho_p - rho_m)

            integrand = HALF * ( upwind*norms(:,1)*unorms(:,1)  +  upwind*norms(:,2)*unorms(:,2)  +  upwind*norms(:,3)*unorms(:,3) )

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irho,integrand)


            !================================
            !       X-MOMENTUM FLUX
            !================================
            upwind = -wave*(rhou_p - rhou_m)

            integrand = HALF * ( upwind*norms(:,1)*unorms(:,1)  +  upwind*norms(:,2)*unorms(:,2)  +  upwind*norms(:,3)*unorms(:,3) )

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhou,integrand)


            !================================
            !       Y-MOMENTUM FLUX
            !================================
            upwind = -wave*(rhov_p - rhov_m)

            integrand = HALF * ( upwind*norms(:,1)*unorms(:,1)  +  upwind*norms(:,2)*unorms(:,2)  +  upwind*norms(:,3)*unorms(:,3) )

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhov,integrand)

            !================================
            !       Z-MOMENTUM FLUX
            !================================
            upwind = -wave*(rhow_p - rhow_m)

            integrand = HALF * ( upwind*norms(:,1)*unorms(:,1)  +  upwind*norms(:,2)*unorms(:,2)  +  upwind*norms(:,3)*unorms(:,3) )

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhow,integrand)

            !================================
            !          ENERGY FLUX
            !================================
            upwind = -wave*(rhoE_p - rhoE_m)

            integrand = HALF * ( upwind*norms(:,1)*unorms(:,1)  +  upwind*norms(:,2)*unorms(:,2)  +  upwind*norms(:,3)*unorms(:,3) )

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhoE,integrand)

        end associate

    end subroutine compute
    !***************************************************************************************************************************************













end module LINEULER_LaxFriedrichs_flux_real
