module bc_euler_wall
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: TWO, HALF, ZERO

    use type_bc,            only: bc_t
    use type_solverdata,    only: solverdata_t
    use type_mesh,          only: mesh_t
    use type_properties,    only: properties_t
    use type_seed,          only: seed_t

    use mod_DNAD_tools,     only: compute_seed
    use mod_integrate,      only: integrate_boundary_scalar_flux
    use mod_interpolate,    only: interpolate_face
    use DNAD_D
    implicit none
    


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_t) :: euler_wall_t

    contains
        procedure :: compute    !> bc implementation
    end type euler_wall_t
    !-------------------------------------------------------------------------------------------




contains

    !> Specialized compute routine for Euler Slip-Wall Boundary Condition
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
        class(euler_wall_t),            intent(inout)   :: self
        type(mesh_t),                   intent(in)      :: mesh(:)
        type(solverdata_t),             intent(inout)   :: sdata
        class(properties_t),            intent(inout)   :: prop
        integer(ik),                    intent(in)      :: idom
        integer(ik),                    intent(in)      :: ielem
        integer(ik),                    intent(in)      :: iface
        integer(ik),                    intent(in)      :: iblk

        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE

        integer(ik)     :: iface_p, ineighbor, idonor
        type(seed_t)    :: seed

        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(idom)%faces(ielem,iface)%gq%face%nnodes)   ::  &
                        rho_m,  rhou_m, rhov_m, rhow_m, rhoE_m, p_m, flux, flux_x, flux_y, flux_z,  &
                        rhou_bc, rhov_bc, rhow_bc, rhoE_bc, u_bc, v_bc, w_bc, u_m, v_m, w_m, p_bc

        real(rk)    :: gam_m


        idonor = 0


        !
        ! Get equation indices
        !
        irho  = prop%get_eqn_index("rho")
        irhou = prop%get_eqn_index("rhou")
        irhov = prop%get_eqn_index("rhov")
        irhow = prop%get_eqn_index("rhow")
        irhoE = prop%get_eqn_index("rhoE")

        !
        ! Get seed element for derivatives
        !
        seed = compute_seed(mesh,idom,ielem,iface,idonor,iblk)



        associate (norms => mesh(idom)%faces(ielem,iface)%norm, unorms => mesh(idom)%faces(ielem,iface)%unorm, faces => mesh(idom)%faces, q => sdata%q)



            !
            ! Interpolate interior solution to quadrature nodes
            !
            call interpolate_face(mesh,q,idom,ielem,iface,irho, rho_m, seed)
            call interpolate_face(mesh,q,idom,ielem,iface,irhou,rhou_m,seed)
            call interpolate_face(mesh,q,idom,ielem,iface,irhov,rhov_m,seed)
            call interpolate_face(mesh,q,idom,ielem,iface,irhow,rhow_m,seed)
            call interpolate_face(mesh,q,idom,ielem,iface,irhoE,rhoE_m,seed)


            !
            ! Compute interior pressure
            !
            call prop%fluid%compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,p_m)
            p_bc = p_m



            !
            ! Initialize arrays
            !
            flux_x = p_bc
            flux_y = p_bc
            flux_z = p_bc
            flux_x = ZERO
            flux_y = ZERO
            flux_z = ZERO



            !
            ! Mass Flux
            !
            flux_x = ZERO
            flux_y = ZERO
            flux_z = ZERO
            flux = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh(idom)%faces(ielem,iface),sdata,idom,irho,iblk,flux)


            !
            ! Add pressure flux to momentum equation
            !
            flux_x = p_bc
            flux_y = ZERO
            flux_z = ZERO
            flux = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh(idom)%faces(ielem,iface),sdata,idom,irhou,iblk,flux)



            !
            ! Add pressure flux to momentum equation
            !
            flux_x = ZERO
            flux_y = p_bc
            flux_z = ZERO

            flux = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh(idom)%faces(ielem,iface),sdata,idom,irhov,iblk,flux)



            !
            ! Add pressure flux to momentum equation
            !
            flux_x = ZERO
            flux_y = ZERO
            flux_z = p_bc

            flux = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh(idom)%faces(ielem,iface),sdata,idom,irhow,iblk,flux)


            !
            ! Energy Flux
            !
            flux_x = ZERO
            flux_y = ZERO
            flux_z = ZERO

            flux = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh(idom)%faces(ielem,iface),sdata,idom,irhoE,iblk,flux)

        end associate

    end subroutine






end module bc_euler_wall
