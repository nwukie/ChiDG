module EULER_Roe_flux
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF,PI, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,ZERO, &
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

    use EULER_properties,       only: EULER_properties_t
    implicit none

    private



    
    !> Implementation of Roe's approximate Riemann solver.
    !!
    !! The formulation used here is from the reference:
    !!   J. Blazek,"Computational Fluid Dynamics: Principles and Applications"
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !------------------------------------------------------------------------------
    type, extends(boundary_flux_t), public :: EULER_Roe_flux_t

    contains

        procedure  :: compute

    end type EULER_Roe_flux_t
    !*******************************************************************************










contains






    !> Compute Roe approximate Riemann upwind flux
    !! 
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!---------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,face_info,function_info)
        class(EULER_Roe_flux_t),            intent(in)      :: self
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


        integer(ik)     :: idom,  ielem,  iface
        integer(ik)     :: ifcn,  idonor, iblk

        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(face_info%idomain_l)%faces(face_info%ielement_l,face_info%iface)%gq%face%nnodes)    ::   &
                        rho_m,      rho_p,                                          &
                        rhou_m,     rhou_p,                                         &
                        rhov_m,     rhov_p,                                         &
                        rhow_m,     rhow_p,                                         &
                        rhoe_m,     rhoe_p,                                         &
                        rhor_p,     rhot_p,                                         &
                        p_m,        p_p,                                            &
                        un_m,       un_p,                                           &
                        a_m,        a_p,                                            &
                        gam_m,      gam_p,                                          &
                        H_m,        H_p,                                            &
                        rtil, util, vtil, wtil, vmagtil, Htil, ctil, qtil2,         &
                        integrand,  upwind,     wave,                               &
                        C1,  C2_a, C2_b,  C3,                                       &
                        u_m, v_m, w_m,                                              &
                        u_p, v_p, w_p,                                              &
                        vmag_p, vmag_m,                                             &
                        delr,   delp,   delvmag, delu, delv, delw,                  &
                        sqrt_rhom, sqrt_rhop, sqrt_rhom_plus_rhop, ctil2, invrho_m, invrho_p


        real(rk), dimension(mesh(face_info%idomain_l)%faces(face_info%ielement_l,face_info%iface)%gq%face%nnodes)    :: &
                        x_m, y_m, z_m, r_m, theta_m, theta_p

        real(rk)    :: theta_offset



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
        idonor = function_info%idonor
        iblk   = function_info%iblk



        associate (norms => mesh(idom)%faces(ielem,iface)%norm, unorms=> mesh(idom)%faces(ielem,iface)%unorm, faces => mesh(idom)%faces, q => sdata%q)

            !
            ! Interpolate solution to quadrature nodes
            !
            call interpolate_face(mesh,face_info,q, irho,  rho_m, LOCAL)
            call interpolate_face(mesh,face_info,q, irho,  rho_p, NEIGHBOR)

            call interpolate_face(mesh,face_info,q, irhou, rhou_m, LOCAL)
            call interpolate_face(mesh,face_info,q, irhou, rhou_p, NEIGHBOR)

            call interpolate_face(mesh,face_info,q, irhov, rhov_m, LOCAL)
            call interpolate_face(mesh,face_info,q, irhov, rhov_p, NEIGHBOR)

            call interpolate_face(mesh,face_info,q, irhow, rhow_m, LOCAL)
            call interpolate_face(mesh,face_info,q, irhow, rhow_p, NEIGHBOR)

            call interpolate_face(mesh,face_info,q, irhoE, rhoE_m, LOCAL)
            call interpolate_face(mesh,face_info,q, irhoE, rhoE_p, NEIGHBOR)



!            idom_m  = face%idomain_l
!            ielem_m = face%ielement_l
!            iface_m = face%iface
!
!
!            if ( mesh(idom_m)%faces(ielem_m, iface_m)%periodic_type == 'cylindrical' ) then
!                theta_offset = (PI/180._rk)*mesh(idom_m)%faces(ielem_m,iface_m)%chimera_offset_theta
!
!                x_m     = mesh(idom_m)%faces(ielem_m,iface_m)%quad_pts(:)%c1_
!                y_m     = mesh(idom_m)%faces(ielem_m,iface_m)%quad_pts(:)%c2_
!                r_m     = sqrt(x_m**TWO + y_m**TWO)
!                theta_m = atan2(y_m,x_m)
!
!
!                theta_p = theta_m + theta_offset
!
!                !
!                ! Transform 'p' cartesian vector to 'p' cylindrical at theta_p
!                !
!                ! rhor and rhotheta are constant from face to face
!                !
!                rhor_p =  cos(theta_p)*rhou_p + sin(theta_p)*rhov_p ! rhor
!                rhot_p = -sin(theta_p)*rhou_p + cos(theta_p)*rhov_p ! rhotheta
!
!
!                !
!                ! Transform 'p' cylindrical to 'p' cartesian at theta_m
!                !
!                rhou_p = cos(theta_m)*rhor_p - sin(theta_m)*rhot_p
!                rhov_p = sin(theta_m)*rhor_p + cos(theta_m)*rhot_p
!
!            end if
!








            !
            ! Compute pressure and gamma
            !
            call prop%fluid%compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,p_m)
            call prop%fluid%compute_pressure(rho_p,rhou_p,rhov_p,rhow_p,rhoE_p,p_p)
            call prop%fluid%compute_gamma(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,gam_m)
            call prop%fluid%compute_gamma(rho_p,rhou_p,rhov_p,rhow_p,rhoE_p,gam_p)

            invrho_m = ONE/rho_m
            invrho_p = ONE/rho_p

            !
            ! Compute enthalpy
            !
            !H_m = (rhoE_m + p_m)/rho_m
            !H_p = (rhoE_p + p_p)/rho_p
            H_m = (rhoE_m + p_m)*invrho_m
            H_p = (rhoE_p + p_p)*invrho_p


            !
            ! Compute velocity components
            !
            !u_m = rhou_m/rho_m
            !v_m = rhov_m/rho_m
            !w_m = rhow_m/rho_m
            u_m = rhou_m*invrho_m
            v_m = rhov_m*invrho_m
            w_m = rhow_m*invrho_m
            vmag_m = u_m*unorms(:,1) + v_m*unorms(:,2) + w_m*unorms(:,3)

            u_p = rhou_p*invrho_p
            v_p = rhov_p*invrho_p
            w_p = rhow_p*invrho_p
            vmag_p = u_p*(unorms(:,1)) + v_p*(unorms(:,2)) + w_p*(unorms(:,3))


            !
            ! Compute Roe-averaged variables
            !
            sqrt_rhom = sqrt(rho_m)
            sqrt_rhop = sqrt(rho_p)
            sqrt_rhom_plus_rhop = sqrt_rhom + sqrt_rhop
            rtil =  sqrt(rho_p * rho_m)                                       ! Roe-averaged density
            util = (sqrt_rhom*u_m + sqrt_rhop*u_p) / (sqrt_rhom_plus_rhop)    ! Roe-averaged u-velocity
            vtil = (sqrt_rhom*v_m + sqrt_rhop*v_p) / (sqrt_rhom_plus_rhop)    ! Roe-averaged v-velocity
            wtil = (sqrt_rhom*w_m + sqrt_rhop*w_p) / (sqrt_rhom_plus_rhop)    ! Roe-averaged w-velocity
            Htil = (sqrt_rhom*H_m + sqrt_rhop*H_p) / (sqrt_rhom_plus_rhop)    ! Roe-averaged Enthalpy

            vmagtil = util*unorms(:,1) + vtil*unorms(:,2) + wtil*unorms(:,3)  ! Magnitude of Roe-averaged velocity in the face normal direction
            qtil2   = util**TWO + vtil**TWO + wtil**TWO

            !& HARDCODED GAMMA
            ctil = sqrt((1.4_rk - ONE)*(Htil - HALF*qtil2))                   ! Roe-averaged speed of sound
            ctil2 = ctil**TWO



            !
            ! Compute jump terms
            !
            delr    = (rho_m - rho_p)
            delu    = (u_m - u_p)
            delv    = (v_m - v_p)
            delw    = (w_m - w_p)
            delvmag = (vmag_m - vmag_p)
            delp    = (p_m - p_p)


            C1   = abs(vmagtil - ctil)*( delp - rtil*ctil*delvmag)/(TWO*(ctil2))
            C2_a = abs(vmagtil)*(delr - delp/(ctil2))
            C2_b = abs(vmagtil)*rtil
            C3   = abs(vmagtil + ctil)*( delp + rtil*ctil*delvmag)/(TWO*(ctil2))


            integrand = delr
            integrand = ZERO


            !================================
            !       MASS FLUX
            !================================
            upwind = C1 + C2_a + C3

            integrand = HALF*(upwind*norms(:,1)*unorms(:,1) + upwind*norms(:,2)*unorms(:,2) + upwind*norms(:,3)*unorms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irho,integrand)


            !================================
            !       X-MOMENTUM FLUX
            !================================
            upwind = C1*(util - ctil*unorms(:,1))  +  C2_a*util  +  C2_b*(delu - delvmag*unorms(:,1))  +  C3*(util + ctil*unorms(:,1))

            integrand = HALF*(upwind*norms(:,1)*unorms(:,1) + upwind*norms(:,2)*unorms(:,2) + upwind*norms(:,3)*unorms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhou,integrand)


            !================================
            !       Y-MOMENTUM FLUX
            !================================
            upwind = C1*(vtil - ctil*unorms(:,2))  +  C2_a*vtil  +  C2_b*(delv - delvmag*unorms(:,2))  +  C3*(vtil + ctil*unorms(:,2))

            integrand = HALF*(upwind*norms(:,1)*unorms(:,1) + upwind*norms(:,2)*unorms(:,2) + upwind*norms(:,3)*unorms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhov,integrand)

            !================================
            !       Z-MOMENTUM FLUX
            !================================
            upwind = C1*(wtil - ctil*unorms(:,3))  +  C2_a*wtil  +  C2_b*(delw - delvmag*unorms(:,3))  +  C3*(wtil + ctil*unorms(:,3))

            integrand = HALF*(upwind*norms(:,1)*unorms(:,1) + upwind*norms(:,2)*unorms(:,2) + upwind*norms(:,3)*unorms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhow,integrand)

            !================================
            !          ENERGY FLUX
            !================================
            upwind = C1*(Htil - ctil*vmagtil)  +  C2_a*(qtil2/TWO)  +  C2_b*(util*delu + vtil*delv + wtil*delw - vmagtil*delvmag)  +  C3*(Htil + ctil*vmagtil)

            integrand = HALF*(upwind*norms(:,1)*unorms(:,1) + upwind*norms(:,2)*unorms(:,2) + upwind*norms(:,3)*unorms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhoE,integrand)

        end associate

    end subroutine compute
    !**********************************************************************************************













end module EULER_Roe_flux
