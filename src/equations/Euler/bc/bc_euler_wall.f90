module bc_euler_wall
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: TWO, HALF, ZERO, ME

    use type_bc,            only: bc_t
    use type_solverdata,    only: solverdata_t
    use type_mesh,          only: mesh_t
    use type_properties,    only: properties_t
    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t

    use mod_integrate,      only: integrate_boundary_scalar_flux
    use mod_interpolate,    only: interpolate
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
    !*******************************************************************************************




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
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,face,fcn)
        class(euler_wall_t),            intent(inout)   :: self
        type(mesh_t),                   intent(in)      :: mesh(:)
        type(solverdata_t),             intent(inout)   :: sdata
        class(properties_t),            intent(inout)   :: prop
        type(face_info_t),              intent(in)      :: face
        type(function_info_t),          intent(in)      :: fcn


        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE

        integer(ik)             :: idom, ielem, iface, idonor

        ! Storage at quadrature nodes
        !type(AD_D), dimension(mesh(face%idomain_l)%faces(face%ielement_l,face%iface)%gq%face%nnodes)   ::  &
        type(AD_D), allocatable, dimension(:)   ::  &
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




        associate (norms => mesh(idom)%faces(ielem,iface)%norm, unorms => mesh(idom)%faces(ielem,iface)%unorm, faces => mesh(idom)%faces, q => sdata%q)



            !
            ! Interpolate interior solution to quadrature nodes
            !
            rho_m  = interpolate(mesh,sdata,face,fcn,irho,  'value', ME)
            rhou_m = interpolate(mesh,sdata,face,fcn,irhou, 'value', ME)
            rhov_m = interpolate(mesh,sdata,face,fcn,irhov, 'value', ME)
            rhow_m = interpolate(mesh,sdata,face,fcn,irhow, 'value', ME)
            rhoE_m = interpolate(mesh,sdata,face,fcn,irhoE, 'value', ME)



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

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irho,integrand)


            !
            ! Add pressure flux to momentum equation
            !
            flux_x = p_bc
            flux_y = ZERO
            flux_z = ZERO
            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irhou,integrand)



            !
            ! Add pressure flux to momentum equation
            !
            flux_x = ZERO
            flux_y = p_bc
            flux_z = ZERO

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irhov,integrand)



            !
            ! Add pressure flux to momentum equation
            !
            flux_x = ZERO
            flux_y = ZERO
            flux_z = p_bc

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irhow,integrand)


            !
            ! Energy Flux
            !
            flux_x = ZERO
            flux_y = ZERO
            flux_z = ZERO

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irhoE,integrand)

        end associate

    end subroutine compute
    !*****************************************************************************************************






end module bc_euler_wall
