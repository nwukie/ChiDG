module LINEULER_LaxFriedrichs_flux_imag
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF,FOUR,ZERO, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX, &
                                      LOCAL, NEIGHBOR

    use atype_boundary_flux,    only: boundary_flux_t
    use type_mesh,              only: mesh_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    use type_seed,              only: seed_t
    use type_face_info,         only: face_info_t
    use type_function_info,     only: function_info_t

    use mod_interpolate,        only: interpolate_face
    use mod_integrate,          only: integrate_boundary_scalar_flux
    use mod_DNAD_tools
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
    !-------------------------------------------------------------------------------------------
    type, extends(boundary_flux_t), public :: LINEULER_LaxFriedrichs_flux_imag_t

    contains

        procedure  :: compute

    end type LINEULER_LaxFriedrichs_flux_imag_t
    !********************************************************************************************










contains






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/16/2017
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,face_info,function_info)
        class(LINEULER_LaxFriedrichs_flux_imag_t),  intent(in)      :: self
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

        integer(ik)     :: idom, ielem, iface
        integer(ik)     :: ifcn, idonor, iblk, igq


        real(rk)    :: a_c


        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(face_info%idomain)%faces(face_info%ielement,face_info%iface)%gq%face%nnodes)    :: &
                        rho_m,      rho_p,                                          &
                        rhou_m,     rhou_p,                                         &
                        rhov_m,     rhov_p,                                         &
                        rhow_m,     rhow_p,                                         &
                        rhoe_m,     rhoe_p,                                         &
                        integrand,  upwind,                                         &
                        wave 

        real(rk), dimension(mesh(face_info%idomain)%faces(face_info%ielement,face_info%iface)%gq%face%nnodes)    :: &
                        un, wave_c


        !===========================================================================
        ! NOTE: var_m signifies "minus" and would indicate a local element variable
        !       var_p signifies "plus"  and would indicate a neighbor element variable
        !===========================================================================
        irho  = prop%get_eqn_index("rho_i")
        irhou = prop%get_eqn_index("rhou_i")
        irhov = prop%get_eqn_index("rhov_i")
        irhow = prop%get_eqn_index("rhow_i")
        irhoE = prop%get_eqn_index("rhoE_i")

        idom  = face_info%idomain
        ielem = face_info%ielement
        iface = face_info%iface
        
        ifcn   = function_info%ifcn
        idonor = function_info%idonor
        iblk   = function_info%iblk


        associate (norms => mesh(idom)%faces(ielem,iface)%norm, unorms=> mesh(idom)%faces(ielem,iface)%unorm, faces => mesh(idom)%faces, q => sdata%q)

            !
            ! Interpolate solution to quadrature nodes
            !
            call interpolate_face(mesh,face_info,q, irho,  rho_m,  LOCAL)
            call interpolate_face(mesh,face_info,q, irho,  rho_p,  NEIGHBOR)

            call interpolate_face(mesh,face_info,q, irhou, rhou_m, LOCAL)
            call interpolate_face(mesh,face_info,q, irhou, rhou_p, NEIGHBOR)

            call interpolate_face(mesh,face_info,q, irhov, rhov_m, LOCAL)
            call interpolate_face(mesh,face_info,q, irhov, rhov_p, NEIGHBOR)

            call interpolate_face(mesh,face_info,q, irhow, rhow_m, LOCAL)
            call interpolate_face(mesh,face_info,q, irhow, rhow_p, NEIGHBOR)

            call interpolate_face(mesh,face_info,q, irhoE, rhoE_m, LOCAL)
            call interpolate_face(mesh,face_info,q, irhoE, rhoE_p, NEIGHBOR)




            wave = rho_m
            do igq = 1,size(wave)
                wave(igq)%x_ad_  = ZERO
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

            integrand = HALF*(upwind*norms(:,1)*unorms(:,1)  +  upwind*norms(:,2)*unorms(:,2)  +  upwind*norms(:,3)*unorms(:,3) )

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irho,integrand)


            !================================
            !       X-MOMENTUM FLUX
            !================================
            upwind = -wave*(rhou_p - rhou_m)

            integrand = HALF*(upwind*norms(:,1)*unorms(:,1)  +  upwind*norms(:,2)*unorms(:,2)  +  upwind*norms(:,3)*unorms(:,3) )

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhou,integrand)


            !================================
            !       Y-MOMENTUM FLUX
            !================================
            upwind = -wave*(rhov_p - rhov_m)

            integrand = HALF*(upwind*norms(:,1)*unorms(:,1)  +  upwind*norms(:,2)*unorms(:,2)  +  upwind*norms(:,3)*unorms(:,3) )

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhov,integrand)

            !================================
            !       Z-MOMENTUM FLUX
            !================================
            upwind = -wave*(rhow_p - rhow_m)

            integrand = HALF*(upwind*norms(:,1)*unorms(:,1)  +  upwind*norms(:,2)*unorms(:,2)  +  upwind*norms(:,3)*unorms(:,3) )

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhow,integrand)

            !================================
            !          ENERGY FLUX
            !================================
            upwind = -wave*(rhoE_p - rhoE_m)

            integrand = HALF*(upwind*norms(:,1)*unorms(:,1)  +  upwind*norms(:,2)*unorms(:,2)  +  upwind*norms(:,3)*unorms(:,3) )

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhoE,integrand)

        end associate

    end subroutine compute
    !**********************************************************************************************************************













end module LINEULER_LaxFriedrichs_flux_imag
