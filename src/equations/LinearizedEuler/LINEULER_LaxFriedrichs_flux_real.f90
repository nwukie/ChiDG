module LINEULER_LaxFriedrichs_flux_real
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF,ZERO, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX, &
                                      LOCAL, NEIGHBOR

    use atype_boundary_flux,    only: boundary_flux_t
    use type_mesh,              only: mesh_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    use type_seed,              only: seed_t
    use type_face_location,     only: face_location_t

    use mod_interpolate,        only: interpolate_face
    use mod_integrate,          only: integrate_boundary_scalar_flux
    use mod_DNAD_tools
    use DNAD_D

    use LINEULER_properties,    only: LINEULER_properties_t
    use mod_linearized_euler
    implicit none

    private

    type, extends(boundary_flux_t), public :: LINEULER_LaxFriedrichs_flux_real_t

    contains
        procedure  :: compute
    end type LINEULER_LaxFriedrichs_flux_real_t










contains






    !==========================================================
    !
    !   Boundary Flux routine for Euler
    !
    !===========================================================
    subroutine compute(self,mesh,sdata,prop,idom,ielem,iface,iblk,idonor)
        class(LINEULER_LaxFriedrichs_flux_real_t),  intent(in)      :: self
        type(mesh_t),                       intent(in)      :: mesh(:)
        type(solverdata_t),                 intent(inout)   :: sdata
        class(properties_t),                intent(inout)   :: prop
        integer(ik),                        intent(in)      :: idom, ielem, iface, iblk
        integer(ik),                        intent(in)      :: idonor

        ! Equation indices
        integer(ik)     :: irho
        integer(ik)     :: irhou
        integer(ik)     :: irhov
        integer(ik)     :: irhow
        integer(ik)     :: irhoe

        integer(ik)             :: igq
        type(seed_t)            :: seed
        type(face_location_t)   :: face

        real(rk)        :: gam_m, gam_p

        real(rk)    :: umag, a_c, wave_c


        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(idom)%faces(ielem,iface)%gq%face%nnodes)    :: &
                        rho_m,      rho_p,                                        &
                        rhou_m,     rhou_p,                                       &
                        rhov_m,     rhov_p,                                       &
                        rhow_m,     rhow_p,                                       &
                        rhoe_m,     rhoe_p,                                       &
                        flux,       upwind,     wave, test_a, test_b


        !===========================================================================
        ! NOTE: var_m signifies "minus" and would indicate a local element variable
        !       var_p signifies "plus"  and would indicate a neighbor element variable
        !===========================================================================
        irho  = prop%get_eqn_index("rho_r")
        irhou = prop%get_eqn_index("rhou_r")
        irhov = prop%get_eqn_index("rhov_r")
        irhow = prop%get_eqn_index("rhow_r")
        irhoE = prop%get_eqn_index("rhoE_r")


        face%idomain  = idom
        face%ielement = ielem
        face%iface    = iface



        !
        ! Compute element for linearization
        !
        seed = compute_seed(mesh,idom,ielem,iface,idonor,iblk)

        associate (norms => mesh(idom)%faces(ielem,iface)%norm, unorms=> mesh(idom)%faces(ielem,iface)%unorm, faces => mesh(idom)%faces, q => sdata%q)

            !
            ! Interpolate solution to quadrature nodes
            !
            call interpolate_face(mesh,sdata%q,idom,ielem,iface, irho,  rho_m,  seed, LOCAL)
            call interpolate_face(mesh,sdata%q,idom,ielem,iface, irho,  rho_p,  seed, NEIGHBOR)

            call interpolate_face(mesh,sdata%q,idom,ielem,iface, irhou, rhou_m, seed, LOCAL)
            call interpolate_face(mesh,sdata%q,idom,ielem,iface, irhou, rhou_p, seed, NEIGHBOR)

            call interpolate_face(mesh,sdata%q,idom,ielem,iface, irhov, rhov_m, seed, LOCAL)
            call interpolate_face(mesh,sdata%q,idom,ielem,iface, irhov, rhov_p, seed, NEIGHBOR)

            call interpolate_face(mesh,sdata%q,idom,ielem,iface, irhow, rhow_m, seed, LOCAL)
            call interpolate_face(mesh,sdata%q,idom,ielem,iface, irhow, rhow_p, seed, NEIGHBOR)

            call interpolate_face(mesh,sdata%q,idom,ielem,iface, irhoE, rhoE_m, seed, LOCAL)
            call interpolate_face(mesh,sdata%q,idom,ielem,iface, irhoE, rhoE_p, seed, NEIGHBOR)


            wave = rho_m
            do igq = 1,size(wave)
                wave(igq)%x_ad_ = ZERO
                wave(igq)%xp_ad_ = ZERO
            end do

            !--------------------------------------
            !  Compute wave speeds
            !--------------------------------------

            ! Compute speed of sound
            a_c = sqrt(abs(gam * pbar / rho_c))



            ! Compute wave speeds
            wave_c = abs(umag) + a_c

            wave = wave_c


            !================================
            !       MASS FLUX
            !================================
            upwind = wave*(rho_m - rho_p)

            flux = HALF*(upwind*norms(:,1)*unorms(:,1)  +  upwind*norms(:,2)*unorms(:,2)  +  upwind*norms(:,3)*unorms(:,3) )

            call integrate_boundary_scalar_flux(mesh,sdata,face,irho,iblk,idonor,seed,flux)


            !================================
            !       X-MOMENTUM FLUX
            !================================
            upwind = wave*(rhou_m - rhou_p)

            flux = HALF*( upwind*norms(:,1)*unorms(:,1)  +  upwind*norms(:,2)*unorms(:,2)  +  upwind*norms(:,3)*unorms(:,3) )

            call integrate_boundary_scalar_flux(mesh,sdata,face,irhou,iblk,idonor,seed,flux)


            !================================
            !       Y-MOMENTUM FLUX
            !================================
            upwind = wave*(rhov_m - rhov_p)

            flux = HALF*( upwind*norms(:,1)*unorms(:,1)  +  upwind*norms(:,2)*unorms(:,2)  +  upwind*norms(:,3)*unorms(:,3) )

            call integrate_boundary_scalar_flux(mesh,sdata,face,irhov,iblk,idonor,seed,flux)

            !================================
            !       Z-MOMENTUM FLUX
            !================================
            upwind = wave*(rhow_m - rhow_p)

            flux = HALF*( upwind*norms(:,1)*unorms(:,1) +  upwind*norms(:,2)*unorms(:,2)  +  upwind*norms(:,3)*unorms(:,3) )

            call integrate_boundary_scalar_flux(mesh,sdata,face,irhow,iblk,idonor,seed,flux)

            !================================
            !          ENERGY FLUX
            !================================
            upwind = wave*(rhoE_m - rhoE_p)

            flux = HALF*( upwind*norms(:,1)*unorms(:,1) + upwind*norms(:,2)*unorms(:,2) + upwind*norms(:,3)*unorms(:,3) )

            call integrate_boundary_scalar_flux(mesh,sdata,face,irhoE,iblk,idonor,seed,flux)

        end associate

    end subroutine













end module LINEULER_LaxFriedrichs_flux_real
