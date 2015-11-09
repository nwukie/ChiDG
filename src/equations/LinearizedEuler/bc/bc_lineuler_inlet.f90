module bc_lineuler_inlet
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ONE, TWO, HALF, ZERO, LOCAL
    use type_bc,            only: bc_t
    use type_solverdata,    only: solverdata_t
    use type_mesh,          only: mesh_t
    use type_properties,    only: properties_t
    use type_seed,          only: seed_t
    use type_face_location, only: face_location_t

    use mod_DNAD_tools,     only: compute_seed
    use mod_integrate,      only: integrate_boundary_scalar_flux
    use mod_interpolate,    only: interpolate_face
    use DNAD_D
    
    use LINEULER_properties,   only: LINEULER_properties_t
    use mod_linearized_euler
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_t) :: lineuler_inlet_t

    contains
        procedure :: compute    !> bc implementation
    end type lineuler_inlet_t
    !-------------------------------------------------------------------------------------------




contains

    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[in]      ielem   Index of the element being computed
    !!  @param[in]      iface   Index of the face being computed
    !!  @param[in]      iblk    Index of the linearization block being computed
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,idom,ielem,iface,iblk)
        class(lineuler_inlet_t),     intent(inout)   :: self
        type(mesh_t),                   intent(in)      :: mesh(:)
        type(solverdata_t),             intent(inout)   :: sdata
        class(properties_t),            intent(inout)   :: prop
        integer(ik),                    intent(in)      :: idom
        integer(ik),                    intent(in)      :: ielem
        integer(ik),                    intent(in)      :: iface
        integer(ik),                    intent(in)      :: iblk

        ! Equation indices
        integer(ik)     :: irho_r, irhou_r, irhov_r, irhow_r, irhoE_r
        integer(ik)     :: irho_i, irhou_i, irhov_i, irhow_i, irhoE_i

        integer(ik)             :: idonor
        type(seed_t)            :: seed
        type(face_location_t)   :: face

        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(idom)%faces(ielem,iface)%gq%face%nnodes)   ::  &
                        rho_r,  rhou_r, rhov_r, rhow_r, rhoE_r,        &
                        rho_i,  rhou_i, rhov_i, rhow_i, rhoE_i,        &
                        flux_x, flux_y, flux_z, flux



        idonor = 0


        !
        ! Get equation indices
        !
        irho_r  = prop%get_eqn_index("rho_r")
        irhou_r = prop%get_eqn_index("rhou_r")
        irhov_r = prop%get_eqn_index("rhov_r")
        irhow_r = prop%get_eqn_index("rhow_r")
        irhoE_r = prop%get_eqn_index("rhoE_r")


        irho_i  = prop%get_eqn_index("rho_i")
        irhou_i = prop%get_eqn_index("rhou_i")
        irhov_i = prop%get_eqn_index("rhov_i")
        irhow_i = prop%get_eqn_index("rhow_i")
        irhoE_i = prop%get_eqn_index("rhoE_i")




        face%idomain  = idom
        face%ielement = ielem
        face%iface    = iface



        !
        ! Get seed element for derivatives
        !
        seed = compute_seed(mesh,idom,ielem,iface,idonor,iblk)


        associate (norms => mesh(idom)%faces(ielem,iface)%norm, unorms => mesh(idom)%faces(ielem,iface)%unorm, faces => mesh(idom)%faces, q => sdata%q)

            !
            ! Interpolate interior solution to quadrature nodes
            !
            call interpolate_face(mesh,q,idom,ielem,iface,irho_r, rho_r, seed, LOCAL)
            call interpolate_face(mesh,q,idom,ielem,iface,irhou_r,rhou_r,seed, LOCAL)
            call interpolate_face(mesh,q,idom,ielem,iface,irhov_r,rhov_r,seed, LOCAL)
            call interpolate_face(mesh,q,idom,ielem,iface,irhow_r,rhow_r,seed, LOCAL)
            call interpolate_face(mesh,q,idom,ielem,iface,irhoE_r,rhoE_r,seed, LOCAL)


            call interpolate_face(mesh,q,idom,ielem,iface,irho_i, rho_i, seed, LOCAL)
            call interpolate_face(mesh,q,idom,ielem,iface,irhou_i,rhou_i,seed, LOCAL)
            call interpolate_face(mesh,q,idom,ielem,iface,irhov_i,rhov_i,seed, LOCAL)
            call interpolate_face(mesh,q,idom,ielem,iface,irhow_i,rhow_i,seed, LOCAL)
            call interpolate_face(mesh,q,idom,ielem,iface,irhoE_i,rhoE_i,seed, LOCAL)


            rho_r  = ZERO
            rhou_r = ZERO
            rhov_r = ZERO
            rhow_r = ZERO
            rhoE_r = ZERO


            rho_i  = ZERO
            rhou_i = ZERO
            rhov_i = ZERO
            rhow_i = ZERO
            rhoE_i = ZERO


            rho_r  = -4.2602822704736e-6_rk
            rhou_r = -5.2125010830968e-9_rk
            rhov_r = ZERO
            rhow_r = ZERO
            rhoE_r = 1.25_rk






            !=================================================
            ! Mass flux
            !=================================================
            flux_x = rho_x_rho  * rho_r  + &
                     rho_x_rhou * rhou_r + &
                     rho_x_rhov * rhov_r + &
                     rho_x_rhoE * rhoE_r
            flux_y = rho_y_rho  * rho_r  + &
                     rho_y_rhou * rhou_r + &
                     rho_y_rhov * rhov_r + &
                     rho_y_rhoE * rhoE_r
            flux_z = rhow_r

            flux = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,irho_r,iblk,idonor,seed,flux)





            flux_x = rho_x_rho  * rho_i  + &
                     rho_x_rhou * rhou_i + &
                     rho_x_rhov * rhov_i + &
                     rho_x_rhoE * rhoE_i
            flux_y = rho_y_rho  * rho_i  + &
                     rho_y_rhou * rhou_i + &
                     rho_y_rhov * rhov_i + &
                     rho_y_rhoE * rhoE_i
            flux_z = rhow_r

            flux = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,irho_i,iblk,idonor,seed,flux)









            !=================================================
            ! x-momentum flux
            !=================================================
            flux_x = rhou_x_rho  * rho_r  + &
                     rhou_x_rhou * rhou_r + &
                     rhou_x_rhov * rhov_r + &
                     rhou_x_rhoE * rhoE_r
            flux_y = rhou_y_rho  * rho_r  + &
                     rhou_y_rhou * rhou_r + &
                     rhou_y_rhov * rhov_r + &
                     rhou_y_rhoE * rhoE_r
            flux_z = ZERO

            flux = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,irhou_r,iblk,idonor,seed,flux)


            flux_x = rhou_x_rho  * rho_i  + &
                     rhou_x_rhou * rhou_i + &
                     rhou_x_rhov * rhov_i + &
                     rhou_x_rhoE * rhoE_i
            flux_y = rhou_y_rho  * rho_i  + &
                     rhou_y_rhou * rhou_i + &
                     rhou_y_rhov * rhov_i + &
                     rhou_y_rhoE * rhoE_i
            flux_z = ZERO

            flux = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,irhou_i,iblk,idonor,seed,flux)







            !=================================================
            ! y-momentum flux
            !=================================================

            flux_x = rhov_x_rho  * rho_r  + &
                     rhov_x_rhou * rhou_r + &
                     rhov_x_rhov * rhov_r + &
                     rhov_x_rhoE * rhoE_r
            flux_y = rhov_y_rho  * rho_r  + &
                     rhov_y_rhou * rhou_r + &
                     rhov_y_rhov * rhov_r + &
                     rhov_y_rhoE * rhoE_r
            flux_z = ZERO

            flux = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,irhov_r,iblk,idonor,seed,flux)


            flux_x = rhov_x_rho  * rho_i  + &
                     rhov_x_rhou * rhou_i + &
                     rhov_x_rhov * rhov_i + &
                     rhov_x_rhoE * rhoE_i
            flux_y = rhov_y_rho  * rho_i  + &
                     rhov_y_rhou * rhou_i + &
                     rhov_y_rhov * rhov_i + &
                     rhov_y_rhoE * rhoE_i
            flux_z = ZERO

            flux = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,irhov_i,iblk,idonor,seed,flux)











!            !=================================================
!            ! z-momentum flux
!            !=================================================
!            flux_x = (rho_m * w_m * u_m)
!            flux_y = (rho_m * w_m * v_m)
!            flux_z = (rho_m * w_m * w_m) + p_m
!
!            flux = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)
!
!            !call integrate_boundary_scalar_flux(mesh(idom)%faces(ielem,iface),sdata,idom,irhow,iblk,flux)
!            call integrate_boundary_scalar_flux(mesh,sdata,face,irhow,iblk,idonor,seed,flux)
!

            !=================================================
            ! Energy flux
            !=================================================


            flux_x = rhoE_x_rho  * rho_r  + &
                     rhoE_x_rhou * rhou_r + &
                     rhoE_x_rhov * rhov_r + &
                     rhoE_x_rhoE * rhoE_r
            flux_y = rhoE_y_rho  * rho_r  + &
                     rhoE_y_rhou * rhou_r + &
                     rhoE_y_rhov * rhov_r + &
                     rhoE_y_rhoE * rhoE_r
            flux_z = ZERO

            flux = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,irhoE_r,iblk,idonor,seed,flux)


            flux_x = rhoE_x_rho  * rho_i  + &
                     rhoE_x_rhou * rhou_i + &
                     rhoE_x_rhov * rhov_i + &
                     rhoE_x_rhoE * rhoE_i
            flux_y = rhoE_y_rho  * rho_i  + &
                     rhoE_y_rhou * rhou_i + &
                     rhoE_y_rhov * rhov_i + &
                     rhoE_y_rhoE * rhoE_i
            flux_z = ZERO

            flux = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,irhoE_i,iblk,idonor,seed,flux)











        end associate

    end subroutine






end module bc_lineuler_inlet
