module bc_euler_pressureoutlet
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ONE, TWO, HALF, LOCAL
    use type_bc,            only: bc_t
    use type_solverdata,    only: solverdata_t
    use type_mesh,          only: mesh_t
    use type_properties,    only: properties_t
    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t

    use mod_integrate,      only: integrate_boundary_scalar_flux
    use mod_interpolate,    only: interpolate_face
    use DNAD_D
    
    use EULER_properties,   only: EULER_properties_t
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !----------------------------------------------------------------------------------------
    type, public, extends(bc_t) :: euler_pressureoutlet_t

    contains

        procedure   :: init_spec    !< user-specialized boundary condition initialization
        procedure   :: compute      !< boundary condition function implementation

    end type euler_pressureoutlet_t
    !****************************************************************************************




contains

    !> Specialized initialization for the boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine init_spec()



    end subroutine init_spec
    !*******************************************************************************************







    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[in]      ielem   Index of the element being computed
    !!  @param[in]      iface   Index of the face being computed
    !!  @param[in]      iblk    Index of the linearization block being computed
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,face,flux)
        class(euler_pressureoutlet_t),  intent(inout)   :: self
        type(mesh_t),                   intent(in)      :: mesh(:)
        type(solverdata_t),             intent(inout)   :: sdata
        class(properties_t),            intent(inout)   :: prop
        type(face_info_t),              intent(in)      :: face
        type(function_info_t),          intent(in)      :: flux


!        integer(ik),                    intent(in)      :: idom
!        integer(ik),                    intent(in)      :: ielem
!        integer(ik),                    intent(in)      :: iface
!        integer(ik),                    intent(in)      :: iblk

        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE

!        type(face_indices_t)    :: face
        integer(ik)     :: idom, ielem, iface, idonor, iblk

        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(face%idomain)%faces(face%ielement,face%iface)%gq%face%nnodes)   ::  &
                        rho_m,  rhou_m, rhov_m, rhow_m, rhoE_m,             &
                        flux_x, flux_y, flux_z, integrand,                  &
                        u_m,    v_m,    w_m,                                &
                        H_bc,   rhoE_bc

        real(rk)    :: gam_m

        real(rk),   dimension(mesh(face%idomain)%faces(face%ielement,face%iface)%gq%face%nnodes)   :: p_bc


        idonor = 0

        !
        ! Get equation indices
        !
        irho  = prop%get_eqn_index("rho")
        irhou = prop%get_eqn_index("rhou")
        irhov = prop%get_eqn_index("rhov")
        irhow = prop%get_eqn_index("rhow")
        irhoE = prop%get_eqn_index("rhoE")



        idom  = face%idomain
        ielem = face%ielement
        iface = face%iface

        iblk   = flux%iblk



        !
        ! Set back pressure
        !
        p_bc = 93000._rk
        !p_bc = 107000._rk



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
            ! Compute velocity components
            !
            u_m = rhou_m/rho_m
            v_m = rhov_m/rho_m
            w_m = rhow_m/rho_m


            !& HARDCODED GAMMA
            gam_m = 1.4_rk


            !
            ! Compute boundary condition energy and enthalpy
            !
            rhoE_bc = p_bc/(gam_m - ONE) + (rho_m/TWO)*(u_m*u_m + v_m*v_m + w_m*w_m)
            H_bc = (rhoE_bc + p_bc)/rho_m

            !=================================================
            ! Mass flux
            !=================================================
            flux_x = (rho_m * u_m)
            flux_y = (rho_m * v_m)
            flux_z = (rho_m * w_m)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irho,integrand)

            !=================================================
            ! x-momentum flux
            !=================================================
            flux_x = (rho_m * u_m * u_m) + p_bc
            flux_y = (rho_m * u_m * v_m)
            flux_z = (rho_m * u_m * w_m)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irhou,integrand)

            !=================================================
            ! y-momentum flux
            !=================================================
            flux_x = (rho_m * v_m * u_m)
            flux_y = (rho_m * v_m * v_m) + p_bc
            flux_z = (rho_m * v_m * w_m)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irhov,integrand)

            !=================================================
            ! z-momentum flux
            !=================================================
            flux_x = (rho_m * w_m * u_m)
            flux_y = (rho_m * w_m * v_m)
            flux_z = (rho_m * w_m * w_m) + p_bc

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irhow,integrand)


            !=================================================
            ! Energy flux
            !=================================================
            flux_x = (rho_m * u_m * H_bc)
            flux_y = (rho_m * v_m * H_bc)
            flux_z = (rho_m * w_m * H_bc)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irhoE,integrand)


        end associate

    end subroutine compute
    !**********************************************************************************************






end module bc_euler_pressureoutlet
