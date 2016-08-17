module EULER_LaxFriedrichs_flux
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX, &
                                      ME, NEIGHBOR

    use type_boundary_flux,     only: boundary_flux_t
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

    use EULER_properties,       only: EULER_properties_t
    implicit none

    private


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!
    !----------------------------------------------------------------------------------
    type, extends(boundary_flux_t), public :: EULER_LaxFriedrichs_flux_t

    contains

        procedure  :: compute

    end type EULER_LaxFriedrichs_flux_t
    !**********************************************************************************










contains




    !>  Compute Lax-Friedrichs upwind flux
    !!
    !!  Dissipation = -alpha(u_m - u_p)
    !!
    !!  Alpha is the maximum wave speed
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/16/2016
    !!
    !!------------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,face_info,function_info)
        class(EULER_LaxFriedrichs_flux_t),  intent(in)      :: self
        type(mesh_t),                       intent(in)      :: mesh(:)
        type(solverdata_t),                 intent(inout)   :: sdata
        class(properties_t),                intent(inout)   :: prop
        type(face_info_t),                  intent(in)      :: face_info
        type(function_info_t),              intent(in)      :: function_info

        ! Equation indices
        integer(ik)     :: irho
        integer(ik)     :: irhou
        integer(ik)     :: irhov
        integer(ik)     :: irhow
        integer(ik)     :: irhoe

        integer(ik)     :: idom, ielem,  iface
        integer(ik)     :: ifcn, idonor, iblk


        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(face_info%idomain_l)%faces(face_info%ielement_l,face_info%iface)%gq%face%nnodes)    :: &
                        rho_m,      rho_p,                                        &
                        rhou_m,     rhou_p,                                       &
                        rhov_m,     rhov_p,                                       &
                        rhow_m,     rhow_p,                                       &
                        rhoe_m,     rhoe_p,                                       &
                        p_m,        p_p,                                          &
                        un_m,       un_p,                                         &
                        a_m,        a_p,                                          &
                        wave_m,     wave_p,                                       &
                        upwind,     wave,                                         &
                        gam_m,      gam_p,                                        &
                        integrand

        real(rk), dimension(mesh(face_info%idomain_l)%faces(face_info%ielement_l,face_info%iface)%gq%face%nnodes)    :: &
                        norm_mag

        !===========================================================================
        ! NOTE: var_m signifies "minus" and would indicate a local element variable
        !       var_p signifies "plus"  and would indicate a neighbor element variable
        !===========================================================================
        irho  = prop%get_eqn_index("rho")
        irhou = prop%get_eqn_index("rhou")
        irhov = prop%get_eqn_index("rhov")
        irhow = prop%get_eqn_index("rhow")
        irhoE = prop%get_eqn_index("rhoE")

        idom  = face_info%idomain_l
        ielem = face_info%ielement_l
        iface = face_info%iface

        ifcn   = function_info%ifcn
        idonor = function_info%idepend
        iblk   = function_info%iblk



        associate (norms => mesh(idom)%faces(ielem,iface)%norm, unorms=> mesh(idom)%faces(ielem,iface)%unorm, faces => mesh(idom)%faces, q => sdata%q)

            !
            ! Interpolate solution to quadrature nodes
            !
            call interpolate_face(mesh,face_info,function_info,q, irho,  rho_m, 'value',  ME)
            call interpolate_face(mesh,face_info,function_info,q, irho,  rho_p, 'value',  NEIGHBOR)

            call interpolate_face(mesh,face_info,function_info,q, irhou, rhou_m, 'value', ME)
            call interpolate_face(mesh,face_info,function_info,q, irhou, rhou_p, 'value', NEIGHBOR)

            call interpolate_face(mesh,face_info,function_info,q, irhov, rhov_m, 'value', ME)
            call interpolate_face(mesh,face_info,function_info,q, irhov, rhov_p, 'value', NEIGHBOR)

            call interpolate_face(mesh,face_info,function_info,q, irhow, rhow_m, 'value', ME)
            call interpolate_face(mesh,face_info,function_info,q, irhow, rhow_p, 'value', NEIGHBOR)

            call interpolate_face(mesh,face_info,function_info,q, irhoE, rhoE_m, 'value', ME)
            call interpolate_face(mesh,face_info,function_info,q, irhoE, rhoE_p, 'value', NEIGHBOR)


            !
            ! Compute pressure and gamma
            !
            call prop%fluid%compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,p_m)
            call prop%fluid%compute_pressure(rho_p,rhou_p,rhov_p,rhow_p,rhoE_p,p_p)
            call prop%fluid%compute_gamma(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,gam_m)
            call prop%fluid%compute_gamma(rho_p,rhou_p,rhov_p,rhow_p,rhoE_p,gam_p)


            !
            ! Compute normal velocities: dot-product vector projection along unit-normal direction
            !
            un_m = unorms(:,1)*(rhou_m/rho_m) + unorms(:,2)*(rhov_m/rho_m) + unorms(:,3)*(rhow_m/rho_m)
            un_p = unorms(:,1)*(rhou_p/rho_p) + unorms(:,2)*(rhov_p/rho_p) + unorms(:,3)*(rhow_p/rho_p)

            
            !
            ! Compute speed of sound
            !
            a_m = sqrt(abs(gam_m * p_m / rho_m))
            a_p = sqrt(abs(gam_p * p_p / rho_p))


            !
            ! Compute wave speeds
            !
            wave_m = abs(un_m) + a_m
            wave_p = abs(un_p) + a_p
            wave   = max(wave_m,wave_p)


            norm_mag = sqrt(norms(:,1)**TWO + norms(:,2)**TWO + norms(:,3)**TWO)

            !================================
            !       MASS FLUX
            !================================
            upwind = -wave*(rho_p - rho_m)

            integrand = HALF * upwind * norm_mag

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irho,integrand)


            !================================
            !       X-MOMENTUM FLUX
            !================================
            upwind = -wave*(rhou_p - rhou_m)

            integrand = HALF * upwind * norm_mag

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhou,integrand)


            !================================
            !       Y-MOMENTUM FLUX
            !================================
            upwind = -wave*(rhov_p - rhov_m)

            integrand = HALF * upwind * norm_mag

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhov,integrand)

            !================================
            !       Z-MOMENTUM FLUX
            !================================
            upwind = -wave*(rhow_p - rhow_m)

            integrand = HALF * upwind * norm_mag

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhow,integrand)

            !================================
            !          ENERGY FLUX
            !================================
            upwind = -wave*(rhoE_p - rhoE_m)

            integrand = HALF * upwind * norm_mag

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhoE,integrand)

        end associate

    end subroutine compute
    !*******************************************************************************************













end module EULER_LaxFriedrichs_flux
