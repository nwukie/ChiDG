module bc_euler_giles_inlet
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ONE, TWO, FOUR, HALF, ZERO, LOCAL

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
    !!  @date   2/16/2016
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_t) :: euler_giles_inlet_t

    contains

        procedure   :: add_options          !< Add boundary condition options.
        procedure   :: init_spec            !< Specialized bc initialization.
        procedure   :: compute              !< bc function implementation.

    end type euler_giles_inlet_t
    !*******************************************************************************************




contains




    !>  Add options for total pressure/temperature boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/16/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine add_options(self)
        class(euler_giles_inlet_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('euler_giles_inlet')


        !
        ! Add functions
        !
        call self%bcproperties%add('TotalPressure',   'Required')
        call self%bcproperties%add('TotalTemperature','Required')
        call self%bcproperties%add('nx',              'Required')
        call self%bcproperties%add('ny',              'Required')
        call self%bcproperties%add('nz',              'Required')
        call self%bcproperties%add('periodicity',     'Required')


        !
        ! Set default total pressure/temperature
        !
        call self%set_fcn_option('TotalPressure',    'val', 110000.0)
        call self%set_fcn_option('TotalTemperature', 'val', 300.0)

        !
        ! Set default angle
        !
        call self%set_fcn_option('nx', 'val', 1._rk)
        call self%set_fcn_option('ny', 'val', 0._rk)
        call self%set_fcn_option('nz', 'val', 0._rk)


    end subroutine add_options
    !********************************************************************************************










    !>  Specialized initialization for the boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/16/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init_spec(self,mesh,iface)
        class(euler_giles_inlet_t), intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh(:)
        integer(ik),                intent(in)      :: iface



    end subroutine init_spec
    !*********************************************************************************************
        






















    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/16/2016
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
        class(euler_giles_inlet_t), intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh(:)
        type(solverdata_t),         intent(inout)   :: sdata
        class(properties_t),        intent(inout)   :: prop
        type(face_info_t),          intent(in)      :: face
        type(function_info_t),      intent(in)      :: flux


        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE


        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(face%idomain)%faces(face%ielement,face%iface)%gq%face%nnodes)   ::  &
                        rho_m,  rhou_m, rhov_m, rhow_m, rhoE_m, p_m,        &
                        flux_x, flux_y, flux_z, integrand,                  &
                        u_m,    v_m,    w_m,                                &
                        u_b,    v_b,    w_b,                                &
                        t_b,    p_b,    rho_b, rhoE_b,                      &
                        vmag2_m, vmag, H_b, Ht, Rplus, c_i, gam_m, a, b, c, cb_plus, cb_minus, c_b, vmag_b, M_b

        real(rk), dimension(mesh(face%idomain)%faces(face%ielement,face%iface)%gq%face%nnodes) :: TT, PT, nx, ny, nz, periodicity

        integer(ik)     :: iface_p, ineighbor, idonor
        integer(ik)     :: idom, ielem, iface, iblk

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
        iblk  = flux%iblk


        associate (norms => mesh(idom)%faces(ielem,iface)%norm, unorms => mesh(idom)%faces(ielem,iface)%unorm, faces => mesh(idom)%faces, &
                coords => mesh(idom)%faces(ielem,iface)%quad_pts,    q => sdata%q,      time => sdata%t )


            !
            ! Get boundary condition Total Temperature, Total Pressure, and normal vector
            !
            PT = self%bcproperties%compute("TotalPressure",        time, coords)
            TT = self%bcproperties%compute("TotalTemperature",     time, coords)
            nx = self%bcproperties%compute("nx",                   time, coords)
            ny = self%bcproperties%compute("ny",                   time, coords)
            nz = self%bcproperties%compute("nz",                   time, coords)
            periodicity = self%bcproperties%compute("periodicity", time, coords)



            !
            ! Get points across boundary at current element span.
            !
            call compute_dft_points(mesh,idom,ielem,iface,periodicity)




            !
            ! Interpolate primitive variables at specified points across boundary. To be DFT'd.
            !
            call interpolate_boundary(mesh,face,q,irho ,rho_p ,points)
            call interpolate_boundary(mesh,face,q,irhou,rhou_p,points)
            call interpolate_boundary(mesh,face,q,irhov,rhov_p,points)
            call interpolate_boundary(mesh,face,q,irhow,rhow_p,points)
            call interpolate_boundary(mesh,face,q,irhoE,rhoE_p,points)




            !
            ! Compute span DFT of primitive variables
            !
            rho_modes  = dft(rho_p )
            rhou_modes = dft(rhou_p)
            rhov_modes = dft(rhov_p)
            rhow_modes = dft(rhow_p)
            rhoE_modes = dft(rhoE_p)








            !
            ! Interpolate interior solution to quadrature nodes
            !
            call interpolate_face(mesh,face,q,irho, rho_m, LOCAL)
            call interpolate_face(mesh,face,q,irhou,rhou_m,LOCAL)
            call interpolate_face(mesh,face,q,irhov,rhov_m,LOCAL)
            call interpolate_face(mesh,face,q,irhow,rhow_m,LOCAL)
            call interpolate_face(mesh,face,q,irhoE,rhoE_m,LOCAL)

            call prop%fluid%compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,p_m)
            call prop%fluid%compute_gamma(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,gam_m)





            !
            ! Compute velocity components
            !
            u_m = rhou_m/rho_m
            v_m = rhov_m/rho_m
            w_m = rhow_m/rho_m


            !
            ! Compute interior speed of sound
            !
            c_i = sqrt(gam_m * p_m / rho_m )


            !
            ! Compute velocity magnitude from interior state
            !
            vmag2_m = (u_m*u_m) + (v_m*v_m) + (w_m*w_m)
            vmag = sqrt(vmag2_m)



            !
            ! Compute interior total enthalpy and R+ characteristic
            !
            Ht    = (p_m / rho_m) * (gam_m/(gam_m - ONE)) + HALF*(vmag2_m)
            Rplus = -vmag - TWO*c_i/(gam_m - ONE)



            !
            ! solve quadratic equation for c_b
            !
            a = ONE + TWO/(gam_m - ONE)
            b = TWO * Rplus
            c = ((gam_m - ONE)/TWO) * (Rplus**TWO - TWO*Ht)
            

            cb_plus  = (-b  +  sqrt(b**TWO - FOUR*a*c))/(TWO*a)
            cb_minus = (-b  -  sqrt(b**TWO - FOUR*a*c))/(TWO*a)

            c_b = max(cb_plus,cb_minus)



            !
            ! Get boundary condition velocity from extrapolated characteristic and computed c_b
            !
            vmag_b = -TWO*c_b/(gam_m - ONE) - Rplus


            !
            ! Compute boundary Mach number
            !
            M_b = vmag_b/c_b



            !
            ! Compute boundary pressure and temperature
            !
            p_b = PT * (ONE + (gam_m - ONE)/TWO * M_b**TWO )**(-gam_m/(gam_m-ONE))
            t_b = TT * (p_b / PT)**((gam_m-ONE)/gam_m)



            !
            ! Compute boundary condition velocity components from imposed direction
            !
            u_b = vmag_b*nx
            v_b = vmag_b*ny
            w_b = vmag_b*nz




            !
            ! Compute boundary condition density from ideal gas law
            !
            select type(prop)
                type is (EULER_properties_t)
                    rho_b = p_b/(t_b*prop%R)
            end select



            !
            ! Compute energy and enthalpy
            !
            rhoE_b = p_b/(gam_m - ONE) + (rho_b/TWO)*( (u_b*u_b) + (v_b*v_b) + (w_b*w_b) )
            H_b    = (rhoE_b + p_b)/rho_b



            !=================================================
            ! Mass flux
            !=================================================
            flux_x = (rho_b * u_b)
            flux_y = (rho_b * v_b)
            flux_z = (rho_b * w_b)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irho,integrand)

            !=================================================
            ! x-momentum flux
            !=================================================
            flux_x = (rho_b * u_b * u_b) + p_b
            flux_y = (rho_b * u_b * v_b)
            flux_z = (rho_b * u_b * w_b)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irhou,integrand)

            !=================================================
            ! y-momentum flux
            !=================================================
            flux_x = (rho_b * v_b * u_b)
            flux_y = (rho_b * v_b * v_b) + p_b
            flux_z = (rho_b * v_b * w_b)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irhov,integrand)

            !=================================================
            ! z-momentum flux
            !=================================================
            flux_x = (rho_b * w_b * u_b)
            flux_y = (rho_b * w_b * v_b)
            flux_z = (rho_b * w_b * w_b) + p_b

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irhow,integrand)


            !=================================================
            ! Energy flux
            !=================================================
            flux_x = (rho_b * u_b * H_b)
            flux_y = (rho_b * v_b * H_b)
            flux_z = (rho_b * w_b * H_b)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,flux,irhoE,integrand)


        end associate

    end subroutine compute
    !*********************************************************************************************






end module bc_euler_giles_inlet
