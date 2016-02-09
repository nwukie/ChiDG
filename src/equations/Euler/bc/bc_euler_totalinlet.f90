module bc_euler_totalinlet
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ONE, TWO, HALF, ZERO, LOCAL

    use type_bc,            only: bc_t
    use type_solverdata,    only: solverdata_t
    use type_mesh,          only: mesh_t
    use type_properties,    only: properties_t
    use type_seed,          only: seed_t
    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t


    use mod_DNAD_tools,     only: compute_seed
    use mod_integrate,      only: integrate_boundary_scalar_flux
    use mod_interpolate,    only: interpolate_face
    use DNAD_D
    
    use EULER_properties,   only: EULER_properties_t
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_t) :: euler_totalinlet_t

    contains

        procedure   :: add_options  !< Add boundary condition options
        procedure   :: compute      !< bc implementation

    end type euler_totalinlet_t
    !-------------------------------------------------------------------------------------------




contains




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine add_options(self)
        class(euler_totalinlet_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('euler_totalinlet')


        !
        ! Add functions
        !
        call self%bcproperties%add('TotalPressure',   'Required')
        call self%bcproperties%add('TotalTemperature','Required')


    end subroutine add_options
    !********************************************************************************************







    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/6/2016
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[in]      ielem   Index of the element being computed
    !!  @param[in]      iface   Index of the face being computed
    !!  @param[in]      iblk    Index of the linearization block being computed
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-------------------------------------------------------------------------------------------
    !subroutine compute(self,mesh,sdata,prop,idom,ielem,iface,iblk)
    subroutine compute(self,mesh,sdata,prop,face,flux)
        class(euler_totalinlet_t),      intent(inout)   :: self
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

        integer(ik)             :: iface_p, ineighbor, idonor
        integer(ik)             :: idom, ielem, iface, iblk
        type(seed_t)            :: seed
!        type(face_indices_t)    :: face

        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(face%idomain)%faces(face%ielement,face%iface)%gq%face%nnodes)   ::  &
                        rho_m,  rhou_m, rhov_m, rhow_m, rhoE_m, p_m,        &
                        flux_x, flux_y, flux_z, integrand,                  &
                        u_m,    v_m,    w_m,                                &
                        u_bc,   v_bc,   w_bc,                               &
                        T_bc,   p_bc,   rho_bc, rhoE_bc,                    &
                        vmag2_m, vmag, H_bc

        real(rk), dimension(mesh(face%idomain)%faces(face%ielement,face%iface)%gq%face%nnodes) :: TT, PT

        !real(rk)    :: gam_m, cp_m, TT, PT, M
        real(rk)    :: gam_m, cp_m, M
        real(rk)    :: norm_bc(3)


        idonor = 0


        !
        ! Get equation indices
        !
        irho  = prop%get_eqn_index("rho")
        irhou = prop%get_eqn_index("rhou")
        irhov = prop%get_eqn_index("rhov")
        irhow = prop%get_eqn_index("rhow")
        irhoE = prop%get_eqn_index("rhoE")



!        face%idomain  = idom
!        face%ielement = ielem
!        face%iface    = iface

        idom  = face%idomain
        ielem = face%ielement
        iface = face%iface

        iblk   = flux%iblk


        !
        ! Get seed element for derivatives
        !
        seed = compute_seed(mesh,idom,ielem,iface,idonor,iblk)


        !associate (norms => mesh(idom)%faces(ielem,iface)%norm, unorms => mesh(idom)%faces(ielem,iface)%unorm, faces => mesh(idom)%faces, q => sdata%q)
        associate (norms => mesh(idom)%faces(ielem,iface)%norm, unorms => mesh(idom)%faces(ielem,iface)%unorm, faces => mesh(idom)%faces, &
                coords => mesh(idom)%faces(ielem,iface)%quad_pts,    q => sdata%q,      time => sdata%t )

        !
        ! Set boundary condition Total Temperature and Total Pressure
        !
        !TT = 300._rk
        !PT = 110000._rk
        PT = self%bcproperties%compute("TotalPressure",   time,coords)
        TT = self%bcproperties%compute("TotalTemperature",time,coords)

        norm_bc = [ONE, ZERO, ZERO]


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

            !
            ! Compute velocity magnitude squared from interior state
            !
            vmag2_m = (u_m*u_m) + (v_m*v_m) + (w_m*w_m)
            vmag = sqrt(vmag2_m)


            !
            ! Compute boundary condition velocity components from imposed direction
            !
            u_bc = vmag*norm_bc(1)
            v_bc = vmag*norm_bc(2)
            w_bc = vmag*norm_bc(3)




            !
            ! Compute boundary condition temperature and pressure
            !
            !& HARDCODED GAMMA. HARDCODED CP
            gam_m = 1.4_rk

            select type(prop)
                type is (EULER_properties_t)
                    cp_m  = (prop%R)*(gam_m/(gam_m-ONE))
            end select

            T_bc = TT - (vmag2_m)/(TWO*cp_m)
            p_bc = PT*((T_bc/TT)**(gam_m/(gam_m-ONE)))

            !gam_m = 1.4_rk
            !cp_m  = 1000._rk
            !M     = 0.4944152_rk
            !T_bc = TT - (gam_m - ONE)*(M**TWO)*(vmag2_m)/(TWO)
            !p_bc = PT*((T_bc/TT)**(gam_m/(gam_m-ONE)))


            !
            ! Compute boundary condition density from ideal gas law
            !
            select type(prop)
                type is (EULER_properties_t)
                    rho_bc = p_bc/(T_bc*prop%R)
            end select




            !
            ! Compute energy and enthalpy
            !
            rhoE_bc = p_bc/(gam_m - ONE) + (rho_bc/TWO)*( (u_bc*u_bc) + (v_bc*v_bc) + (w_bc*w_bc) )
            H_bc    = (rhoE_bc + p_bc)/rho_bc



            !=================================================
            ! Mass flux
            !=================================================
            flux_x = (rho_bc * u_bc)
            flux_y = (rho_bc * v_bc)
            flux_z = (rho_bc * w_bc)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irho,integrand)

            !=================================================
            ! x-momentum flux
            !=================================================
            flux_x = (rho_bc * u_bc * u_bc) + p_bc
            flux_y = (rho_bc * u_bc * v_bc)
            flux_z = (rho_bc * u_bc * w_bc)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irhou,integrand)

            !=================================================
            ! y-momentum flux
            !=================================================
            flux_x = (rho_bc * v_bc * u_bc)
            flux_y = (rho_bc * v_bc * v_bc) + p_bc
            flux_z = (rho_bc * v_bc * w_bc)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irhov,integrand)

            !=================================================
            ! z-momentum flux
            !=================================================
            flux_x = (rho_bc * w_bc * u_bc)
            flux_y = (rho_bc * w_bc * v_bc)
            flux_z = (rho_bc * w_bc * w_bc) + p_bc

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irhow,integrand)


            !=================================================
            ! Energy flux
            !=================================================
            flux_x = (rho_bc * u_bc * H_bc)
            flux_y = (rho_bc * v_bc * H_bc)
            flux_z = (rho_bc * w_bc * H_bc)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irhoE,integrand)


        end associate

    end subroutine






end module bc_euler_totalinlet
