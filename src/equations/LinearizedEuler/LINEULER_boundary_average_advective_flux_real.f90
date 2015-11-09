module LINEULER_boundary_average_advective_flux_real
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES, ZERO, ONE, TWO, HALF, &
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


    integer :: counter = 0

    private



    !> Implementation of the Euler boundary average flux
    !!
    !!  - At a boundary interface, the solution states Q- and Q+ exists on opposite 
    !!    sides of the boundary. The average flux is computed as Favg = 1/2(F(Q-) + F(Q+))
    !!
    !!  @author Nathan A. Wukie
    !!
    !--------------------------------------------------------------------------------
    type, extends(boundary_flux_t), public :: LINEULER_boundary_average_advective_flux_real_t

    contains
        procedure  :: compute

    end type LINEULER_boundary_average_advective_flux_real_t










contains



    !
    !   Boundary Flux routine for Euler
    !
    !----------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,idom,ielem,iface,iblk,idonor)
        class(LINEULER_boundary_average_advective_flux_real_t), intent(in)      :: self
        type(mesh_t),                                   intent(in)      :: mesh(:)
        type(solverdata_t),                             intent(inout)   :: sdata
        class(properties_t),                            intent(inout)   :: prop
        integer(ik),                                    intent(in)      :: idom, ielem, iface, iblk
        integer(ik),                                    intent(in)      :: idonor

        ! Equation indices
        integer(ik)     :: irho
        integer(ik)     :: irhou
        integer(ik)     :: irhov
        integer(ik)     :: irhow
        integer(ik)     :: irhoe



        integer(ik) :: igq




        type(seed_t)            :: seed
        type(face_location_t)   :: face

        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(idom)%faces(ielem,iface)%gq%face%nnodes)    :: &
                        rho_m,      rho_p,                                  &
                        rhou_m,     rhou_p,                                 &
                        rhov_m,     rhov_p,                                 &
                        rhow_m,     rhow_p,                                 &
                        rhoe_m,     rhoe_p,                                 &
                        flux_x_m,   flux_y_m,   flux_z_m,                   &
                        flux_x_p,   flux_y_p,   flux_z_p,                   &
                        flux_x,     flux_y,     flux_z,                     &
                        flux


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




        associate (norms => mesh(idom)%faces(ielem,iface)%norm, faces => mesh(idom)%faces, q => sdata%q)

            !
            ! Compute element for linearization
            !
            seed = compute_seed(mesh,idom,ielem,iface,idonor,iblk)



            !
            ! Interpolate solution to quadrature nodes
            !
            call interpolate_face(mesh,q,idom,ielem,iface, irho,  rho_m,  seed, LOCAL)
            call interpolate_face(mesh,q,idom,ielem,iface, irho,  rho_p,  seed, NEIGHBOR)

            call interpolate_face(mesh,q,idom,ielem,iface, irhou, rhou_m, seed, LOCAL)
            call interpolate_face(mesh,q,idom,ielem,iface, irhou, rhou_p, seed, NEIGHBOR)

            call interpolate_face(mesh,q,idom,ielem,iface, irhov, rhov_m, seed, LOCAL)
            call interpolate_face(mesh,q,idom,ielem,iface, irhov, rhov_p, seed, NEIGHBOR)

            call interpolate_face(mesh,q,idom,ielem,iface, irhow, rhow_m, seed, LOCAL)
            call interpolate_face(mesh,q,idom,ielem,iface, irhow, rhow_p, seed, NEIGHBOR)

            call interpolate_face(mesh,q,idom,ielem,iface, irhoE, rhoE_m, seed, LOCAL)
            call interpolate_face(mesh,q,idom,ielem,iface, irhoE, rhoE_p, seed, NEIGHBOR)








            !================================
            !       MASS FLUX
            !================================
            flux_x_m = rho_x_rho  * rho_m  + &
                       rho_x_rhou * rhou_m + &
                       rho_x_rhov * rhov_m + &
                       rho_x_rhoE * rhoE_m
            flux_y_m = rho_y_rho  * rho_m  + &
                       rho_y_rhou * rhou_m + &
                       rho_y_rhov * rhov_m + &
                       rho_y_rhoE * rhoE_m 
            flux_z_m = rhow_m

            flux_x_p = rho_x_rho  * rho_p  + &
                       rho_x_rhou * rhou_p + &
                       rho_x_rhov * rhov_p + &
                       rho_x_rhoE * rhoE_p
            flux_y_p = rho_y_rho  * rho_p  + &
                       rho_y_rhou * rhou_p + &
                       rho_y_rhov * rhov_p + &
                       rho_y_rhoE * rhoE_p 
            flux_z_p = rhow_p

            flux_x = (flux_x_m + flux_x_p)
            flux_y = (flux_y_m + flux_y_p)
            flux_z = (flux_z_m + flux_z_p)


            ! dot with normal vector
            flux = HALF*(flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face,irho,iblk,idonor,seed,flux)


            !================================
            !       X-MOMENTUM FLUX
            !================================
            flux_x_m = rhou_x_rho  * rho_m  + &
                       rhou_x_rhou * rhou_m + &
                       rhou_x_rhov * rhov_m + &
                       rhou_x_rhoE * rhoE_m
            flux_y_m = rhou_y_rho  * rho_m  + &
                       rhou_y_rhou * rhou_m + &
                       rhou_y_rhov * rhov_m + &
                       rhou_y_rhoE * rhoE_m 
            flux_z_m = ZERO


            flux_x_p = rhou_x_rho  * rho_p  + &
                       rhou_x_rhou * rhou_p + &
                       rhou_x_rhov * rhov_p + &
                       rhou_x_rhoE * rhoE_p
            flux_y_p = rhou_y_rho  * rho_p  + &
                       rhou_y_rhou * rhou_p + &
                       rhou_y_rhov * rhov_p + &
                       rhou_y_rhoE * rhoE_p 
            flux_z_p = ZERO





            flux_x = (flux_x_m + flux_x_p)
            flux_y = (flux_y_m + flux_y_p)
            flux_z = (flux_z_m + flux_z_p)


            ! dot with normal vector
            flux = HALF*(flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face,irhou,iblk,idonor,seed,flux)


            !================================
            !       Y-MOMENTUM FLUX
            !================================
            flux_x_m = rhov_x_rho  * rho_m  + &
                       rhov_x_rhou * rhou_m + &
                       rhov_x_rhov * rhov_m + &
                       rhov_x_rhoE * rhoE_m
            flux_y_m = rhov_y_rho  * rho_m  + &
                       rhov_y_rhou * rhou_m + &
                       rhov_y_rhov * rhov_m + &
                       rhov_y_rhoE * rhoE_m 
            flux_z_m = ZERO



            flux_x_p = rhov_x_rho  * rho_p  + &
                       rhov_x_rhou * rhou_p + &
                       rhov_x_rhov * rhov_p + &
                       rhov_x_rhoE * rhoE_p
            flux_y_p = rhov_y_rho  * rho_p  + &
                       rhov_y_rhou * rhou_p + &
                       rhov_y_rhov * rhov_p + &
                       rhov_y_rhoE * rhoE_p 
            flux_z_p = ZERO





            flux_x = (flux_x_m + flux_x_p)
            flux_y = (flux_y_m + flux_y_p)
            flux_z = (flux_z_m + flux_z_p)


            ! dot with normal vector
            flux = HALF*(flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face,irhov,iblk,idonor,seed,flux)

!            !================================
!            !       Z-MOMENTUM FLUX
!            !================================
!            flux_x_m = (rhow_m*rhou_m)/rho_m
!            flux_y_m = (rhow_m*rhov_m)/rho_m
!            flux_z_m = (rhow_m*rhow_m)/rho_m + p_m
!
!            flux_x_p = (rhow_p*rhou_p)/rho_p
!            flux_y_p = (rhow_p*rhov_p)/rho_p
!            flux_z_p = (rhow_p*rhow_p)/rho_p + p_p
!
!            flux_x = (flux_x_m + flux_x_p)
!            flux_y = (flux_y_m + flux_y_p)
!            flux_z = (flux_z_m + flux_z_p)
!
!
!            ! dot with normal vector
!            flux = HALF*(flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3))
!
!            call integrate_boundary_scalar_flux(mesh,sdata,face,irhow,iblk,idonor,seed,flux)
!
            !================================
            !          ENERGY FLUX
            !================================
            flux_x_m = rhoE_x_rho  * rho_m  + &
                       rhoE_x_rhou * rhou_m + &
                       rhoE_x_rhov * rhov_m + &
                       rhoE_x_rhoE * rhoE_m
            flux_y_m = rhoE_y_rho  * rho_m  + &
                       rhoE_y_rhou * rhou_m + &
                       rhoE_y_rhov * rhov_m + &
                       rhoE_y_rhoE * rhoE_m 
            flux_z_m = ZERO


            flux_x_p = rhoE_x_rho  * rho_p  + &
                       rhoE_x_rhou * rhou_p + &
                       rhoE_x_rhov * rhov_p + &
                       rhoE_x_rhoE * rhoE_p
            flux_y_p = rhoE_y_rho  * rho_p  + &
                       rhoE_y_rhou * rhou_p + &
                       rhoE_y_rhov * rhov_p + &
                       rhoE_y_rhoE * rhoE_p 
            flux_z_p = ZERO





            flux_x = (flux_x_m + flux_x_p)
            flux_y = (flux_y_m + flux_y_p)
            flux_z = (flux_z_m + flux_z_p)


            ! dot with normal vector
            flux = HALF*(flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face,irhoE,iblk,idonor,seed,flux)

        end associate

    end subroutine












end module LINEULER_boundary_average_advective_flux_real
