module EULER_Roe_flux_cacheoptimized
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
    use type_timer,             only: timer_t

    use mod_interpolate,        only: interpolate_face
    use mod_integrate,          only: integrate_boundary_scalar_flux
    use mod_DNAD_tools
    use DNAD_D

    use EULER_properties,       only: EULER_properties_t
    implicit none


    type(timer_t), public   :: EULER_roe_total
    type(timer_t), public   :: EULER_roe_interpolate
    type(timer_t), public   :: EULER_roe_integrate


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
    type, extends(boundary_flux_t), public :: EULER_Roe_flux_cacheoptimized_t

    contains

        procedure  :: compute

    end type EULER_Roe_flux_cacheoptimized_t
    !*******************************************************************************










contains






    !> Compute Roe approximate Riemann upwind flux
    !! 
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!---------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,face_info,function_info)
        class(EULER_Roe_flux_cacheoptimized_t),            intent(in)      :: self
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
        integer(ik)     :: ngq, block, iloop, gq_max, gq_min
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



        call EULER_roe_total%start()



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
            call EULER_roe_interpolate%start()
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
            call EULER_roe_interpolate%stop()



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



            ngq = size(rho_m)
            block = 5


            iloop  = 1
            gq_max = 0
            do while (gq_max < ngq)
                gq_min = 1  +  block*(iloop-1)
                gq_max = min(ngq,block*(iloop))


                !
                ! Compute pressure and gamma
                !
                call prop%fluid%compute_pressure(rho_m(gq_min:gq_max),rhou_m(gq_min:gq_max),rhov_m(gq_min:gq_max),rhow_m(gq_min:gq_max),rhoE_m(gq_min:gq_max),p_m(gq_min:gq_max))
                call prop%fluid%compute_pressure(rho_p(gq_min:gq_max),rhou_p(gq_min:gq_max),rhov_p(gq_min:gq_max),rhow_p(gq_min:gq_max),rhoE_p(gq_min:gq_max),p_p(gq_min:gq_max))
                call prop%fluid%compute_gamma(rho_m(gq_min:gq_max),rhou_m(gq_min:gq_max),rhov_m(gq_min:gq_max),rhow_m(gq_min:gq_max),rhoE_m(gq_min:gq_max),gam_m(gq_min:gq_max))
                call prop%fluid%compute_gamma(rho_p(gq_min:gq_max),rhou_p(gq_min:gq_max),rhov_p(gq_min:gq_max),rhow_p(gq_min:gq_max),rhoE_p(gq_min:gq_max),gam_p(gq_min:gq_max))

                invrho_m(gq_min:gq_max) = ONE/rho_m(gq_min:gq_max)
                invrho_p(gq_min:gq_max) = ONE/rho_p(gq_min:gq_max)

                !
                ! Compute enthalpy
                !
                H_m(gq_min:gq_max) = (rhoE_m(gq_min:gq_max) + p_m(gq_min:gq_max))*invrho_m(gq_min:gq_max)
                H_p(gq_min:gq_max) = (rhoE_p(gq_min:gq_max) + p_p(gq_min:gq_max))*invrho_p(gq_min:gq_max)


                !
                ! Compute velocity components
                !
                u_m(gq_min:gq_max) = rhou_m(gq_min:gq_max)*invrho_m(gq_min:gq_max)
                v_m(gq_min:gq_max) = rhov_m(gq_min:gq_max)*invrho_m(gq_min:gq_max)
                w_m(gq_min:gq_max) = rhow_m(gq_min:gq_max)*invrho_m(gq_min:gq_max)
                vmag_m(gq_min:gq_max) = u_m(gq_min:gq_max)*unorms(gq_min:gq_max,1) + v_m(gq_min:gq_max)*unorms(gq_min:gq_max,2) + w_m(gq_min:gq_max)*unorms(gq_min:gq_max,3)

                u_p(gq_min:gq_max) = rhou_p(gq_min:gq_max)*invrho_p(gq_min:gq_max)
                v_p(gq_min:gq_max) = rhov_p(gq_min:gq_max)*invrho_p(gq_min:gq_max)
                w_p(gq_min:gq_max) = rhow_p(gq_min:gq_max)*invrho_p(gq_min:gq_max)
                vmag_p(gq_min:gq_max) = u_p(gq_min:gq_max)*(unorms(gq_min:gq_max,1)) + v_p(gq_min:gq_max)*(unorms(gq_min:gq_max,2)) + w_p(gq_min:gq_max)*(unorms(gq_min:gq_max,3))


                !
                ! Compute Roe-averaged variables
                !
                sqrt_rhom(gq_min:gq_max) = sqrt(rho_m(gq_min:gq_max))
                sqrt_rhop(gq_min:gq_max) = sqrt(rho_p(gq_min:gq_max))
                sqrt_rhom_plus_rhop(gq_min:gq_max) = sqrt_rhom(gq_min:gq_max) + sqrt_rhop(gq_min:gq_max)
                rtil(gq_min:gq_max) =  sqrt(rho_p(gq_min:gq_max) * rho_m(gq_min:gq_max))                                       ! Roe-averaged density
                util(gq_min:gq_max) = (sqrt_rhom(gq_min:gq_max)*u_m(gq_min:gq_max) + sqrt_rhop(gq_min:gq_max)*u_p(gq_min:gq_max)) / (sqrt_rhom_plus_rhop(gq_min:gq_max))    ! Roe-averaged u-velocity
                vtil(gq_min:gq_max) = (sqrt_rhom(gq_min:gq_max)*v_m(gq_min:gq_max) + sqrt_rhop(gq_min:gq_max)*v_p(gq_min:gq_max)) / (sqrt_rhom_plus_rhop(gq_min:gq_max))    ! Roe-averaged v-velocity
                wtil(gq_min:gq_max) = (sqrt_rhom(gq_min:gq_max)*w_m(gq_min:gq_max) + sqrt_rhop(gq_min:gq_max)*w_p(gq_min:gq_max)) / (sqrt_rhom_plus_rhop(gq_min:gq_max))    ! Roe-averaged w-velocity
                Htil(gq_min:gq_max) = (sqrt_rhom(gq_min:gq_max)*H_m(gq_min:gq_max) + sqrt_rhop(gq_min:gq_max)*H_p(gq_min:gq_max)) / (sqrt_rhom_plus_rhop(gq_min:gq_max))    ! Roe-averaged Enthalpy

                vmagtil(gq_min:gq_max) = util(gq_min:gq_max)*unorms(gq_min:gq_max,1) + vtil(gq_min:gq_max)*unorms(gq_min:gq_max,2) + wtil(gq_min:gq_max)*unorms(gq_min:gq_max,3)  ! Magnitude of Roe-averaged velocity in the face normal direction
                qtil2(gq_min:gq_max)   = util(gq_min:gq_max)**TWO + vtil(gq_min:gq_max)**TWO + wtil(gq_min:gq_max)**TWO

                !& HARDCODED GAMMA
                ctil(gq_min:gq_max) = sqrt((1.4_rk - ONE)*(Htil(gq_min:gq_max) - HALF*qtil2(gq_min:gq_max)))                   ! Roe-averaged speed of sound
                ctil2(gq_min:gq_max) = ctil(gq_min:gq_max)**TWO



                !
                ! Compute jump terms
                !
                delr(gq_min:gq_max)    = (rho_m(gq_min:gq_max) - rho_p(gq_min:gq_max))
                delu(gq_min:gq_max)    = (u_m(gq_min:gq_max) - u_p(gq_min:gq_max))
                delv(gq_min:gq_max)    = (v_m(gq_min:gq_max) - v_p(gq_min:gq_max))
                delw(gq_min:gq_max)    = (w_m(gq_min:gq_max) - w_p(gq_min:gq_max))
                delvmag(gq_min:gq_max) = (vmag_m(gq_min:gq_max) - vmag_p(gq_min:gq_max))
                delp(gq_min:gq_max)    = (p_m(gq_min:gq_max) - p_p(gq_min:gq_max))


                C1(gq_min:gq_max)   = abs(vmagtil(gq_min:gq_max) - ctil(gq_min:gq_max))*( delp(gq_min:gq_max) - rtil(gq_min:gq_max)*ctil(gq_min:gq_max)*delvmag(gq_min:gq_max))/(TWO*(ctil2(gq_min:gq_max)))
                C2_a(gq_min:gq_max) = abs(vmagtil(gq_min:gq_max))*(delr(gq_min:gq_max) - delp(gq_min:gq_max)/(ctil2(gq_min:gq_max)))
                C2_b(gq_min:gq_max) = abs(vmagtil(gq_min:gq_max))*rtil(gq_min:gq_max)
                C3(gq_min:gq_max)   = abs(vmagtil(gq_min:gq_max) + ctil(gq_min:gq_max))*( delp(gq_min:gq_max) + rtil(gq_min:gq_max)*ctil(gq_min:gq_max)*delvmag(gq_min:gq_max))/(TWO*(ctil2(gq_min:gq_max)))


            !integrand = delr
            !integrand = ZERO


            !================================
            !       MASS FLUX
            !================================
            upwind(gq_min:gq_max) = C1(gq_min:gq_max) + C2_a(gq_min:gq_max) + C3(gq_min:gq_max)

            integrand(gq_min:gq_max) = HALF*(upwind(gq_min:gq_max)*norms(gq_min:gq_max,1)*unorms(gq_min:gq_max,1) + upwind(gq_min:gq_max)*norms(gq_min:gq_max,2)*unorms(gq_min:gq_max,2) + upwind(gq_min:gq_max)*norms(gq_min:gq_max,3)*unorms(gq_min:gq_max,3))

                iloop = iloop + 1
            end do

            call EULER_roe_integrate%start()
            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irho,integrand)
            call EULER_roe_integrate%stop()


            !================================
            !       X-MOMENTUM FLUX
            !================================
            iloop  = 1
            gq_max = 0
            do while (gq_max < ngq)
                gq_min = 1  +  block*(iloop-1)
                gq_max = min(ngq,block*(iloop))
            upwind(gq_min:gq_max) = C1(gq_min:gq_max)*(util(gq_min:gq_max) - ctil(gq_min:gq_max)*unorms(gq_min:gq_max,1))  +  C2_a(gq_min:gq_max)*util(gq_min:gq_max)  +  C2_b(gq_min:gq_max)*(delu(gq_min:gq_max) - delvmag(gq_min:gq_max)*unorms(gq_min:gq_max,1))  +  C3(gq_min:gq_max)*(util(gq_min:gq_max) + ctil(gq_min:gq_max)*unorms(gq_min:gq_max,1))

            integrand(gq_min:gq_max) = HALF*(upwind(gq_min:gq_max)*norms(gq_min:gq_max,1)*unorms(gq_min:gq_max,1) + upwind(gq_min:gq_max)*norms(gq_min:gq_max,2)*unorms(gq_min:gq_max,2) + upwind(gq_min:gq_max)*norms(gq_min:gq_max,3)*unorms(gq_min:gq_max,3))

                iloop = iloop + 1
            end do
            call EULER_roe_integrate%start()
            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhou,integrand)
            call EULER_roe_integrate%stop()


            !================================
            !       Y-MOMENTUM FLUX
            !================================
            iloop  = 1
            gq_max = 0
            do while (gq_max < ngq)
                gq_min = 1  +  block*(iloop-1)
                gq_max = min(ngq,block*(iloop))
            upwind(gq_min:gq_max) = C1(gq_min:gq_max)*(vtil(gq_min:gq_max) - ctil(gq_min:gq_max)*unorms(gq_min:gq_max,2))  +  C2_a(gq_min:gq_max)*vtil(gq_min:gq_max)  +  C2_b(gq_min:gq_max)*(delv(gq_min:gq_max) - delvmag(gq_min:gq_max)*unorms(gq_min:gq_max,2))  +  C3(gq_min:gq_max)*(vtil(gq_min:gq_max) + ctil(gq_min:gq_max)*unorms(gq_min:gq_max,2))

            integrand(gq_min:gq_max) = HALF*(upwind(gq_min:gq_max)*norms(gq_min:gq_max,1)*unorms(gq_min:gq_max,1) + upwind(gq_min:gq_max)*norms(gq_min:gq_max,2)*unorms(gq_min:gq_max,2) + upwind(gq_min:gq_max)*norms(gq_min:gq_max,3)*unorms(gq_min:gq_max,3))

                iloop = iloop + 1
            end do
            call EULER_roe_integrate%start()
            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhov,integrand)
            call EULER_roe_integrate%stop()

            !================================
            !       Z-MOMENTUM FLUX
            !================================
            iloop  = 1
            gq_max = 0
            do while (gq_max < ngq)
                gq_min = 1  +  block*(iloop-1)
                gq_max = min(ngq,block*(iloop))
            upwind(gq_min:gq_max) = C1(gq_min:gq_max)*(wtil(gq_min:gq_max) - ctil(gq_min:gq_max)*unorms(gq_min:gq_max,3))  +  C2_a(gq_min:gq_max)*wtil(gq_min:gq_max)  +  C2_b(gq_min:gq_max)*(delw(gq_min:gq_max) - delvmag(gq_min:gq_max)*unorms(gq_min:gq_max,3))  +  C3(gq_min:gq_max)*(wtil(gq_min:gq_max) + ctil(gq_min:gq_max)*unorms(gq_min:gq_max,3))

            integrand(gq_min:gq_max) = HALF*(upwind(gq_min:gq_max)*norms(gq_min:gq_max,1)*unorms(gq_min:gq_max,1) + upwind(gq_min:gq_max)*norms(gq_min:gq_max,2)*unorms(gq_min:gq_max,2) + upwind(gq_min:gq_max)*norms(gq_min:gq_max,3)*unorms(gq_min:gq_max,3))

                iloop = iloop + 1
            end do
            call EULER_roe_integrate%start()
            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhow,integrand)
            call EULER_roe_integrate%stop()

            !================================
            !          ENERGY FLUX
            !================================
            iloop  = 1
            gq_max = 0
            do while (gq_max < ngq)
                gq_min = 1  +  block*(iloop-1)
                gq_max = min(ngq,block*(iloop))
            upwind(gq_min:gq_max) = C1(gq_min:gq_max)*(Htil(gq_min:gq_max) - ctil(gq_min:gq_max)*vmagtil(gq_min:gq_max))  +  C2_a(gq_min:gq_max)*(qtil2(gq_min:gq_max)*HALF)  +  C2_b(gq_min:gq_max)*(util(gq_min:gq_max)*delu(gq_min:gq_max) + vtil(gq_min:gq_max)*delv(gq_min:gq_max) + wtil(gq_min:gq_max)*delw(gq_min:gq_max) - vmagtil(gq_min:gq_max)*delvmag(gq_min:gq_max))  +  C3(gq_min:gq_max)*(Htil(gq_min:gq_max) + ctil(gq_min:gq_max)*vmagtil(gq_min:gq_max))

            integrand(gq_min:gq_max) = HALF*(upwind(gq_min:gq_max)*norms(gq_min:gq_max,1)*unorms(gq_min:gq_max,1) + upwind(gq_min:gq_max)*norms(gq_min:gq_max,2)*unorms(gq_min:gq_max,2) + upwind(gq_min:gq_max)*norms(gq_min:gq_max,3)*unorms(gq_min:gq_max,3))
                iloop = iloop + 1
            end do

            call EULER_roe_integrate%start()
            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irhoE,integrand)
            call EULER_roe_integrate%stop()

        end associate

        call EULER_roe_total%stop()
    end subroutine compute
    !**********************************************************************************************













end module EULER_Roe_flux_cacheoptimized
