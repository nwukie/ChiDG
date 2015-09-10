module eqn_euler
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX
    use atype_equationset,      only: equationset_t
    use type_mesh,              only: mesh_t
    use atype_solverdata,       only: solverdata_t
    use mod_interpolate,        only: interpolate
    use mod_integrate,          only: integrate_volume_flux, integrate_boundary_flux, integrate_boundary_scalar_flux
    use mod_DNAD_tools,         only: compute_neighbor_face, compute_seed_element
    use DNAD_D

    implicit none

    private

    type, extends(equationset_t), public :: euler_e

        ! Equation-set specific data
!        type(fluid_t)   :: fluid

    contains
        ! Must define these procedures in the extended, concrete type
        procedure  :: init
        procedure  :: compute_boundary_average_flux
        procedure  :: compute_boundary_upwind_flux
        procedure  :: compute_volume_flux
        procedure  :: compute_volume_source


    end type euler_e










contains
    !==========================================================
    !
    !   Equation set initialization
    !
    !==========================================================
    subroutine init(self)
        class(euler_e), intent(inout) :: self

        self%neqns   = 5

        ! Allocate equations
        allocate(self%eqns(self%neqns))

        ! Initialize equation parameters
        self%eqns(1)%name = "rho"
        self%eqns(1)%ind  = 1

        self%eqns(2)%name = "rhou"
        self%eqns(2)%ind  = 2

        self%eqns(3)%name = "rhov"
        self%eqns(3)%ind  = 3

        self%eqns(4)%name = "rhow"
        self%eqns(4)%ind  = 4

        self%eqns(5)%name = "rhoE"
        self%eqns(5)%ind  = 5

        ! Initialize equation set parameters
!        self%rgas = 287.058_rk        ! J/(kg*K)

    end subroutine




    !==========================================================
    !
    !   Boundary Flux routine for Euler
    !
    !===========================================================
    subroutine compute_boundary_average_flux(self,mesh,sdata,ielem,iface,iblk)
        class(euler_e),         intent(in)      :: self
        class(mesh_t),          intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iface, iblk

        ! Equation indices
        integer(ik)     :: irho
        integer(ik)     :: irhou
        integer(ik)     :: irhov
        integer(ik)     :: irhow
        integer(ik)     :: irhoe

        integer(ik)     :: iseed, iface_p, ineighbor

        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh%faces(ielem,iface)%gq%face%nnodes)    :: &
                        rho_m,      rho_p,                                  &
                        rhou_m,     rhou_p,                                 &
                        rhov_m,     rhov_p,                                 &
                        rhow_m,     rhow_p,                                 &
                        rhoe_m,     rhoe_p,                                 &
                        p_m,        p_p,                                    &
                        H_m,        H_p,                                    &
                        flux_x_m,   flux_y_m,   flux_z_m,                   &
                        flux_x_p,   flux_y_p,   flux_z_p,                   &
                        flux_x,     flux_y,     flux_z,                     &
                        flux


        !===========================================================================
        ! NOTE: var_m signifies "minus" and would indicate a local element variable
        !       var_p signifies "plus"  and would indicate a neighbor element variable
        !===========================================================================
        irho  = self%get_var("rho")
        irhou = self%get_var("rhou")
        irhov = self%get_var("rhov")
        irhow = self%get_var("rhow")
        irhoE = self%get_var("rhoE")


        associate (norms => mesh%faces(ielem,iface)%norm, faces => mesh%faces, q => sdata%q)

            !> Get neighbor face and seed element for derivatives
            iface_p   = compute_neighbor_face(iface)
            iseed     = compute_seed_element(mesh,ielem,iblk)
            ineighbor = mesh%faces(ielem,iface)%ineighbor




            !> Interpolate solution to quadrature nodes
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
            call compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,p_m)
            call compute_pressure(rho_p,rhou_p,rhov_p,rhow_p,rhoE_p,p_p)

            H_m = (rhoE_m + p_m)/rho_m
            H_p = (rhoE_p + p_p)/rho_p

            !================================
            !       MASS FLUX
            !================================
            flux_x_m = rhou_m
            flux_y_m = rhov_m
            flux_z_m = rhow_m

            flux_x_p = rhou_p
            flux_y_p = rhov_p
            flux_z_p = rhow_p

            flux_x = (flux_x_m + flux_x_p)*HALF*norms(:,1)
            flux_y = (flux_y_m + flux_y_p)*HALF*norms(:,2)
            flux_z = (flux_z_m + flux_z_p)*HALF*norms(:,3)


            call integrate_boundary_flux(mesh%faces(ielem,iface),sdata,irho,iblk,flux_x,flux_y,flux_z)


            !================================
            !       X-MOMENTUM FLUX
            !================================
            flux_x_m = (rhou_m*rhou_m)/rho_m + p_m
            flux_y_m = (rhou_m*rhov_m)/rho_m
            flux_z_m = (rhou_m*rhow_m)/rho_m

            flux_x_p = (rhou_p*rhou_p)/rho_p + p_p
            flux_y_p = (rhou_p*rhov_p)/rho_p
            flux_z_p = (rhou_p*rhow_p)/rho_p

            flux_x = (flux_x_m + flux_x_p)*HALF*norms(:,1)
            flux_y = (flux_y_m + flux_y_p)*HALF*norms(:,2)
            flux_z = (flux_z_m + flux_z_p)*HALF*norms(:,3)

            call integrate_boundary_flux(mesh%faces(ielem,iface),sdata,irhou,iblk,flux_x,flux_y,flux_z)


            !================================
            !       Y-MOMENTUM FLUX
            !================================
            flux_x_m = (rhov_m*rhou_m)/rho_m
            flux_y_m = (rhov_m*rhov_m)/rho_m + p_m
            flux_z_m = (rhov_m*rhow_m)/rho_m

            flux_x_p = (rhov_p*rhou_p)/rho_p
            flux_y_p = (rhov_p*rhov_p)/rho_p + p_p
            flux_z_p = (rhov_p*rhow_p)/rho_p

            flux_x = (flux_x_m + flux_x_p)*HALF*norms(:,1)
            flux_y = (flux_y_m + flux_y_p)*HALF*norms(:,2)
            flux_z = (flux_z_m + flux_z_p)*HALF*norms(:,3)

            call integrate_boundary_flux(mesh%faces(ielem,iface),sdata,irhov,iblk,flux_x,flux_y,flux_z)

            !================================
            !       Z-MOMENTUM FLUX
            !================================
            flux_x_m = (rhow_m*rhou_m)/rho_m
            flux_y_m = (rhow_m*rhov_m)/rho_m
            flux_z_m = (rhow_m*rhow_m)/rho_m + p_m

            flux_x_p = (rhow_p*rhou_p)/rho_p
            flux_y_p = (rhow_p*rhov_p)/rho_p
            flux_z_p = (rhow_p*rhow_p)/rho_p + p_p

            flux_x = (flux_x_m + flux_x_p)*HALF*norms(:,1)
            flux_y = (flux_y_m + flux_y_p)*HALF*norms(:,2)
            flux_z = (flux_z_m + flux_z_p)*HALF*norms(:,3)

            call integrate_boundary_flux(mesh%faces(ielem,iface),sdata,irhow,iblk,flux_x,flux_y,flux_z)

            !================================
            !          ENERGY FLUX
            !================================
            flux_x_m = rhou_m*H_m
            flux_y_m = rhov_m*H_m
            flux_z_m = rhow_m*H_m

            flux_x_p = rhou_p*H_p
            flux_y_p = rhov_p*H_p
            flux_z_p = rhow_p*H_p

            flux_x = (flux_x_m + flux_x_p)*HALF*norms(:,1)
            flux_y = (flux_y_m + flux_y_p)*HALF*norms(:,2)
            flux_z = (flux_z_m + flux_z_p)*HALF*norms(:,3)

            !print*, flux_x%x_ad_
            !stop
            call integrate_boundary_flux(mesh%faces(ielem,iface),sdata,irhoE,iblk,flux_x,flux_y,flux_z)

        end associate

    end subroutine









    !==========================================================
    !
    !   Boundary Flux routine for Euler
    !
    !===========================================================
    subroutine compute_boundary_upwind_flux(self,mesh,sdata,ielem,iface,iblk)
        class(euler_e),         intent(in)      :: self
        class(mesh_t),          intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iface, iblk

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
        irho  = self%get_var("rho")
        irhou = self%get_var("rhou")
        irhov = self%get_var("rhov")
        irhow = self%get_var("rhow")
        irhoE = self%get_var("rhoE")

        !> Get neighbor face and seed element for derivatives
        iface_p   = compute_neighbor_face(iface)
        iseed     = compute_seed_element(mesh,ielem,iblk)
        ineighbor = mesh%faces(ielem,iface)%ineighbor

        associate (norms => mesh%faces(ielem,iface)%norm, unorms=> mesh%faces(ielem,iface)%unorm, faces => mesh%faces, q => sdata%q)

            !> Interpolate solution to quadrature nodes
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
            call compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,p_m)
            call compute_pressure(rho_p,rhou_p,rhov_p,rhow_p,rhoE_p,p_p)

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

            ! Dot product with normal vector
            call integrate_boundary_scalar_flux(mesh%faces(ielem,iface),sdata,irho,iblk,flux)


            !================================
            !       X-MOMENTUM FLUX
            !================================
            upwind = -wave*(rhou_p - rhou_m)

            flux = HALF*(upwind)

            ! Dot product with normal vector
            call integrate_boundary_scalar_flux(mesh%faces(ielem,iface),sdata,irhou,iblk,flux)


            !================================
            !       Y-MOMENTUM FLUX
            !================================
            upwind = -wave*(rhov_p - rhov_m)

            flux = HALF*(upwind)

            ! Dot product with normal vector
            call integrate_boundary_scalar_flux(mesh%faces(ielem,iface),sdata,irhov,iblk,flux)

            !================================
            !       Z-MOMENTUM FLUX
            !================================
            upwind = -wave*(rhow_p - rhow_m)

            flux = HALF*(upwind)

            ! Dot product with normal vector
            call integrate_boundary_scalar_flux(mesh%faces(ielem,iface),sdata,irhow,iblk,flux)

            !================================
            !          ENERGY FLUX
            !================================
            upwind = -wave*(rhoE_p - rhoE_m)

            flux = HALF*(upwind)

            ! Dot product with normal vector
            call integrate_boundary_scalar_flux(mesh%faces(ielem,iface),sdata,irhoE,iblk,flux)

        end associate

    end subroutine














    !==========================================================
    !
    !   Volume Flux routine for Euler
    !
    !===========================================================
    subroutine compute_volume_flux(self,mesh,sdata,ielem,iblk)
        class(euler_e),         intent(in)      :: self
        class(mesh_t),          intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iblk

        ! Equation indices
        !------------------------------------------------------------
        integer(ik)    :: irho
        integer(ik)    :: irhou
        integer(ik)    :: irhov
        integer(ik)    :: irhow
        integer(ik)    :: irhoe


        integer(ik)    :: iseed, i

        type(AD_D), dimension(mesh%elems(ielem)%gq%vol%nnodes)      ::  &
                    rho, rhou, rhov, rhow, rhoE, p, H,                  &
                    flux_x, flux_y, flux_z

        !-------------------------------------------------------------
        irho  = self%get_var("rho")
        irhou = self%get_var("rhou")
        irhov = self%get_var("rhov")
        irhow = self%get_var("rhow")
        irhoE = self%get_var("rhoE")

        !> Get neighbor face and seed element for derivatives
        iseed   = compute_seed_element(mesh,ielem,iblk)


        !> Interpolate solution to quadrature nodes
        call interpolate(mesh%elems,sdata%q,ielem,irho, rho, iseed)
        call interpolate(mesh%elems,sdata%q,ielem,irhou,rhou,iseed)
        call interpolate(mesh%elems,sdata%q,ielem,irhov,rhov,iseed)
        call interpolate(mesh%elems,sdata%q,ielem,irhow,rhow,iseed)
        call interpolate(mesh%elems,sdata%q,ielem,irhoE,rhoE,iseed)


        ! Compute pressure and total enthalpy
        call compute_pressure(rho,rhou,rhov,rhow,rhoE,p)

        H = (rhoE + p)/rho

        !===========================
        !        MASS FLUX
        !===========================
        flux_x = rhou
        flux_y = rhov
        flux_z = rhow

        call integrate_volume_flux(mesh%elems(ielem),sdata,irho,iblk,flux_x,flux_y,flux_z)


        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux_x = (rhou*rhou)/rho  +  p
        flux_y = (rhou*rhov)/rho
        flux_z = (rhou*rhow)/rho

        call integrate_volume_flux(mesh%elems(ielem),sdata,irhou,iblk,flux_x,flux_y,flux_z)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux_x = (rhov*rhou)/rho
        flux_y = (rhov*rhov)/rho  +  p
        flux_z = (rhov*rhow)/rho

        call integrate_volume_flux(mesh%elems(ielem),sdata,irhov,iblk,flux_x,flux_y,flux_z)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        flux_x = (rhow*rhou)/rho
        flux_y = (rhow*rhov)/rho
        flux_z = (rhow*rhow)/rho  +  p

        call integrate_volume_flux(mesh%elems(ielem),sdata,irhow,iblk,flux_x,flux_y,flux_z)

        !============================
        !       ENERGY FLUX
        !============================
        flux_x = rhou*H
        flux_y = rhov*H
        flux_z = rhow*H

        call integrate_volume_flux(mesh%elems(ielem),sdata,irhoE,iblk,flux_x,flux_y,flux_z)

    end subroutine


    !==========================================================
    !
    !   Volume Source routine for Euler
    !
    !===========================================================
    subroutine compute_volume_source(self,mesh,sdata,ielem,iblk)
        class(euler_e),         intent(in)      :: self
        class(mesh_t),          intent(in)      :: mesh
        class(solverdata_t),    intent(inout)   :: sdata
        integer(ik),            intent(in)      :: ielem, iblk




    end subroutine



    !================================================================
    !
    !   Implements the equation of state for computing pressure
    !
    !================================================================
    subroutine compute_pressure(rho,rhou,rhov,rhow,rhoE,p)
        type(AD_D),  intent(in)     :: rho(:), rhou(:), rhov(:), rhow(:), rhoE(:)
        type(AD_D),  intent(inout)  :: p(:)

        real(rk)    :: gam

        gam = 1.4_rk




        p = (gam-ONE)*(rhoE - HALF*rho*((rhou*rhou)/(rho*rho) + (rhov*rhov)/(rho*rho) + (rhow*rhow)/(rho*rho)))




    end subroutine










end module eqn_euler
