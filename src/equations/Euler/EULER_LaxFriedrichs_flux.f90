module EULER_LaxFriedrichs_flux
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX

    use atype_boundary_flux,    only: boundary_flux_t
    use atype_equationset,      only: equationset_t
    use type_mesh,              only: mesh_t
    use atype_solverdata,       only: solverdata_t
    use mod_interpolate,        only: interpolate
    use mod_integrate,          only: integrate_volume_flux, integrate_boundary_flux, integrate_boundary_scalar_flux
    use mod_DNAD_tools,         only: compute_neighbor_face, compute_seed_element
    use DNAD_D

    use type_properties,        only: properties_t
    use EULER_properties,       only: EULER_properties_t
    implicit none

    private

    type, extends(boundary_flux_t), public :: EULER_LaxFriedrichs_flux_t

    contains
        procedure  :: compute
    end type EULER_LaxFriedrichs_flux_t










contains






    !==========================================================
    !
    !   Boundary Flux routine for Euler
    !
    !===========================================================
    subroutine compute(self,mesh,sdata,ielem,iface,iblk,prop)
        class(EULER_LaxFriedrichs_flux_t),  intent(in)      :: self
        class(mesh_t),                      intent(in)      :: mesh
        class(solverdata_t),                intent(inout)   :: sdata
        integer(ik),                        intent(in)      :: ielem, iface, iblk
        class(properties_t),                intent(inout)   :: prop

        ! Equation indices
        integer(ik)     :: irho
        integer(ik)     :: irhou
        integer(ik)     :: irhov
        integer(ik)     :: irhow
        integer(ik)     :: irhoe

        integer(ik)     :: iseed, iface_p, ineighbor

        real(rk)        :: gam_m, gam_p

        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh%faces(ielem,iface)%gq%face%nnodes)    :: &
                        rho_m,      rho_p,                                  &
                        rhou_m,     rhou_p,                                 &
                        rhov_m,     rhov_p,                                 &
                        rhow_m,     rhow_p,                                 &
                        rhoe_m,     rhoe_p,                                 &
                        p_m,        p_p,                                    &
                        un_m,       un_p,                                   &
                        a_m,        a_p,                                    &
                        wave_m,     wave_p,                                 &
                        flux,       upwind,     wave, test_a, test_b


        !===========================================================================
        ! NOTE: var_m signifies "minus" and would indicate a local element variable
        !       var_p signifies "plus"  and would indicate a neighbor element variable
        !===========================================================================
        irho  = prop%get_eqn_index("rho")
        irhou = prop%get_eqn_index("rhou")
        irhov = prop%get_eqn_index("rhov")
        irhow = prop%get_eqn_index("rhow")
        irhoE = prop%get_eqn_index("rhoE")

        ! Get neighbor face and seed element for derivatives
        iface_p   = compute_neighbor_face(iface)
        iseed     = compute_seed_element(mesh,ielem,iblk)
        ineighbor = mesh%faces(ielem,iface)%ineighbor

        associate (norms => mesh%faces(ielem,iface)%norm, unorms=> mesh%faces(ielem,iface)%unorm, faces => mesh%faces, q => sdata%q)

            ! Interpolate solution to quadrature nodes
            call interpolate(faces,q,ielem,    iface,  irho,rho_m,iseed)
            call interpolate(faces,q,ineighbor,iface_p,irho,rho_p,iseed)

            call interpolate(faces,q,ielem,    iface,  irhou,rhou_m,iseed)
            call interpolate(faces,q,ineighbor,iface_p,irhou,rhou_p,iseed)

            call interpolate(faces,q,ielem,    iface,  irhov,rhov_m,iseed)
            call interpolate(faces,q,ineighbor,iface_p,irhov,rhov_p,iseed)

            call interpolate(faces,q,ielem,    iface,  irhow,rhow_m,iseed)
            call interpolate(faces,q,ineighbor,iface_p,irhow,rhow_p,iseed)

            call interpolate(faces,q,ielem,    iface,  irhoE,rhoE_m,iseed)
            call interpolate(faces,q,ineighbor,iface_p,irhoE,rhoE_p,iseed)



            ! Compute pressure and total enthalpy
            call prop%fluid%compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,p_m)
            call prop%fluid%compute_pressure(rho_p,rhou_p,rhov_p,rhow_p,rhoE_p,p_p)

            gam_m = 1.4_rk
            gam_p = 1.4_rk
            !--------------------------------------
            !  Compute wave speeds
            !--------------------------------------
            ! Compute normal velocities: dot-product vector projection along unit-normal direction
            un_m = unorms(:,1)*(rhou_m/rho_m) + unorms(:,2)*(rhov_m/rho_m) + unorms(:,3)*(rhow_m/rho_m)
            un_p = -unorms(:,1)*(rhou_p/rho_p) - unorms(:,2)*(rhov_p/rho_p) - unorms(:,3)*(rhow_p/rho_p)

            ! Compute speed of sound
            a_m = sqrt(abs(gam_m * p_m / rho_m))
            a_p = sqrt(abs(gam_p * p_p / rho_p))



            ! Compute wave speeds
            wave_m = abs(un_m) + a_m
            wave_p = abs(un_p) + a_p
            wave = max(wave_m,wave_p)


            !================================
            !       MASS FLUX
            !================================
            upwind = -wave*(rho_p - rho_m)

            flux = HALF*(upwind)
            !flux = HALF*(upwind)*norms(:,1) + HALF*(upwind)*norms(:,2) + HALF*(upwind)*norms(:,3)

            call integrate_boundary_scalar_flux(mesh%faces(ielem,iface),sdata,irho,iblk,flux)


            !================================
            !       X-MOMENTUM FLUX
            !================================
            upwind = -wave*(rhou_p - rhou_m)

            flux = HALF*(upwind)
            !flux = HALF*(upwind)*norms(:,1) + HALF*(upwind)*norms(:,2) + HALF*(upwind)*norms(:,3)

            call integrate_boundary_scalar_flux(mesh%faces(ielem,iface),sdata,irhou,iblk,flux)


            !================================
            !       Y-MOMENTUM FLUX
            !================================
            upwind = -wave*(rhov_p - rhov_m)

            flux = HALF*(upwind)
            !flux = HALF*(upwind)*norms(:,1) + HALF*(upwind)*norms(:,2) + HALF*(upwind)*norms(:,3)

            call integrate_boundary_scalar_flux(mesh%faces(ielem,iface),sdata,irhov,iblk,flux)

            !================================
            !       Z-MOMENTUM FLUX
            !================================
            upwind = -wave*(rhow_p - rhow_m)

            flux = HALF*(upwind)
            !flux = HALF*(upwind)*norms(:,1) + HALF*(upwind)*norms(:,2) + HALF*(upwind)*norms(:,3)

            call integrate_boundary_scalar_flux(mesh%faces(ielem,iface),sdata,irhow,iblk,flux)

            !================================
            !          ENERGY FLUX
            !================================
            upwind = -wave*(rhoE_p - rhoE_m)

            flux = HALF*(upwind)
            !flux = HALF*(upwind)*norms(:,1) + HALF*(upwind)*norms(:,2) + HALF*(upwind)*norms(:,3)

            call integrate_boundary_scalar_flux(mesh%faces(ielem,iface),sdata,irhoE,iblk,flux)

        end associate

    end subroutine













end module EULER_LaxFriedrichs_flux
