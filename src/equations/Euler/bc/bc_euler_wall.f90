module bc_euler_wall
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: TWO, HALF, ZERO, LOCAL

    use type_bc,            only: bc_t
    use type_solverdata,    only: solverdata_t
    use type_mesh,          only: mesh_t
    use type_properties,    only: properties_t
    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t

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
        procedure   :: add_options
        procedure   :: compute    !> bc implementation
    end type euler_wall_t
    !-------------------------------------------------------------------------------------------




contains



    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_options(self)    
        class(euler_wall_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('euler_wall')


        !
        ! Add functions
        !


        !
        ! Add parameters
        !


    end subroutine add_options
    !******************************************************************************************















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
    subroutine compute(self,mesh,sdata,prop,face,flux)
        class(euler_wall_t),            intent(inout)   :: self
        type(mesh_t),                   intent(in)      :: mesh(:)
        type(solverdata_t),             intent(inout)   :: sdata
        class(properties_t),            intent(inout)   :: prop
        type(face_info_t),              intent(in)      :: face
        type(function_info_t),          intent(in)      :: flux


        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE

        integer(ik)             :: idom, ielem, iface, idonor, iblk

        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(face%idomain_l)%faces(face%ielement_l,face%iface)%gq%face%nnodes)   ::  &
                        rho_m,  rhou_m, rhov_m, rhow_m, rhoE_m, p_m, integrand, flux_x, flux_y, flux_z,  &
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


        idom  = face%idomain_l
        ielem = face%ielement_l
        iface = face%iface

        iblk   = flux%iblk



        associate (norms => mesh(idom)%faces(ielem,iface)%norm, unorms => mesh(idom)%faces(ielem,iface)%unorm, faces => mesh(idom)%faces, q => sdata%q)



            !
            ! Interpolate interior solution to quadrature nodes
            !
            call interpolate_face(mesh,face,q,irho, rho_m, LOCAL)
            call interpolate_face(mesh,face,q,irhou,rhou_m,LOCAL)
            call interpolate_face(mesh,face,q,irhov,rhov_m,LOCAL)
            call interpolate_face(mesh,face,q,irhow,rhow_m,LOCAL)
            call interpolate_face(mesh,face,q,irhoE,rhoE_m,LOCAL)

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
            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irho,integrand)


            !
            ! Add pressure flux to momentum equation
            !
            flux_x = p_bc
            flux_y = ZERO
            flux_z = ZERO
            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irhou,integrand)



            !
            ! Add pressure flux to momentum equation
            !
            flux_x = ZERO
            flux_y = p_bc
            flux_z = ZERO

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irhov,integrand)



            !
            ! Add pressure flux to momentum equation
            !
            flux_x = ZERO
            flux_y = ZERO
            flux_z = p_bc

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irhow,integrand)


            !
            ! Energy Flux
            !
            flux_x = ZERO
            flux_y = ZERO
            flux_z = ZERO

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irhoE,integrand)

        end associate

    end subroutine






end module bc_euler_wall
