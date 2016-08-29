module bc_euler_giles_inlet
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, FOUR, HALF, ZERO

    use type_bc,                only: bc_t
    use type_solverdata,        only: solverdata_t
    use type_mesh,              only: mesh_t
    use type_point,             only: point_t
    use type_properties,        only: properties_t
    use type_face_info,         only: face_info_t
    use type_function_info,     only: function_info_t


    use mod_chidg_interpolate,  only: interpolate
    use mod_integrate,          only: integrate_boundary_scalar_flux
    use mod_dft,                only: dft, idft_mode_points
    use DNAD_D
    
    use EULER_properties,       only: EULER_properties_t
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   2/16/2016
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_t) :: euler_giles_inlet_t

        type(point_t),  allocatable :: dft_points(:)


    contains

        procedure   :: add_options              !< Add boundary condition options.
        procedure   :: init_spec                !< Specialized bc initialization.
        procedure   :: init_boundary_coupling   !< Implement specialized coupling information between elements.
        procedure   :: compute                  !< bc function implementation.

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
        call self%set_fcn_option('TotalPressure',    'val', 110000.0_rk)
        call self%set_fcn_option('TotalTemperature', 'val', 300.0_rk)

        !
        ! Set default angle
        !
        call self%set_fcn_option('nx', 'val', ONE)
        call self%set_fcn_option('ny', 'val', ZERO)
        call self%set_fcn_option('nz', 'val', ZERO)


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
        type(mesh_t),               intent(inout)   :: mesh
        integer(ik),                intent(in)      :: iface

        real(rk)        :: periodicity
        real(rk)        :: zero_time
        type(point_t)   :: zero_point


        !
        ! Get boundary periodicity from bc options.
        !
        zero_time  = ZERO
        call zero_point%set(ZERO,ZERO,ZERO)
        periodicity = self%bcproperties%compute("periodicity", zero_time, zero_point)


        !
        ! Compute dft points
        !
        self%dft_points = compute_dft_points(mesh,self%elems,iface,periodicity)


    end subroutine init_spec
    !*********************************************************************************************
        










    !>  Implement specific boundary coupling.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/16/2016
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------
    subroutine init_boundary_coupling(self,mesh,iface)
        class(euler_giles_inlet_t), intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh
        integer(ik),                intent(in)      :: iface

        integer(ik) :: ielem_bc, ielem_coupled, ielem


        !
        ! Loop through elements. For the current 2D giles, every element on the boundary
        ! is coupled with every other element on the boundary through the Fourier transform.
        !
        do ielem_bc = 1,size(self%elems)


            !
            ! Register all elements as coupled to the current element.
            !
            do ielem_coupled = 1,size(self%elems)

                !
                ! Get block-element index of current ielem_bc.
                !
                ielem = self%elems(ielem_bc) 

                !
                ! Add element index to the coupling for the current element.
                !
                call self%coupled_elems(ielem_bc)%push_back(ielem)

            end do ! ielem_coupled

        end do  !ielem_bc



    end subroutine init_boundary_coupling
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
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,face,fcn)
        class(euler_giles_inlet_t), intent(inout)   :: self
        type(mesh_t),               intent(in)      :: mesh(:)
        type(solverdata_t),         intent(inout)   :: sdata
        class(properties_t),        intent(inout)   :: prop
        type(face_info_t),          intent(in)      :: face
        type(function_info_t),      intent(in)      :: fcn


        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE


        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(face%idomain)%faces(face%ielement,face%iface)%gq%face%nnodes)   ::  &
                        rho_m,  rhou_m, rhov_m, rhow_m, rhoE_m, p_m,        &
                        flux_x, flux_y, flux_z, integrand,                  &
                        u_m,    v_m,    w_m,                                &
                        u_b,    v_b,    w_b,                                &
                        t_b,    p_b,    rho_b, rhoE_b,                      &
                        vmag2_m, vmag, H_b, Ht, Rplus, c_i, gam_m, a, b, c, cb_plus, cb_minus, c_b, vmag_b, M_b, &
                        rho_modes, rhou_modes, rhov_modes, rhow_modes, rhoE_modes

        real(rk), dimension(mesh(face%idomain)%faces(face%ielement,face%iface)%gq%face%nnodes) :: TT, PT, nx, ny, nz, periodicity

        integer(ik)     :: iface_p, ineighbor, idonor
        integer(ik)     :: idom, ielem, iface

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




            !---------------------------------------------------------------------
            !
            !   Decompose boundary into Fourier modes of characteristic variables
            !
            !---------------------------------------------------------------------

            !
            ! Interpolate conservative variables at specified points across boundary. To be DFT'd.
            !
            call interpolate_boundary(mesh,face,fcn,q,irho ,rho_b ,self%dft_points)
            call interpolate_boundary(mesh,face,fcn,q,irhou,rhou_b,self%dft_points)
            call interpolate_boundary(mesh,face,fcn,q,irhov,rhov_b,self%dft_points)
            call interpolate_boundary(mesh,face,fcn,q,irhow,rhow_b,self%dft_points)
            call interpolate_boundary(mesh,face,fcn,q,irhoE,rhoE_b,self%dft_points)

        
            !
            ! Compute primitive variables
            !
            call prop%fluid%compute_pressure(rho_b,rhou_b,rhov_b,rhow_b,rhoE_b,p_b  )
            call prop%fluid%compute_gamma(   rho_b,rhou_b,rhov_b,rhow_b,rhoE_b,gam_b)
            u_b = rhou_b / rho_b
            v_b = rhov_b / rho_b
            w_b = rhow_b / rho_b


            !
            ! Compute speed of sound across boundary
            !
            c_b = sqrt(gam_b * p_b / rho_b )



            !
            ! Compute characteristic variables across boundary
            !
            c1_b = -(c_b**TWO)*rho_b +      ZERO     +       ZERO      +    p_b
            c2_b =         ZERO      +      ZERO     +  rho_b*c_b*v_b  +   ZERO
            c3_b =         ZERO      + rho_b*c_b*u_b +       ZERO      +    p_b
            c4_b =         ZERO      - rho_b*c_b*u_b +       ZERO      +    p_b


            !
            ! Compute DFT of characteristic variables
            !
            call dft(c1_b, c1_real, c1_imag)
            call dft(c2_b, c2_real, c2_imag)
            call dft(c3_b, c3_real, c3_imag)
            call dft(c4_b, c4_real, c4_imag)





            !-------------------------------------------------------
            !
            !   Adjust mean component
            !
            !-------------------------------------------------------


            !
            ! Evaluate mean component of characteristic variables at gq nodes
            !
            imode = 1
            call idft_mode_points(ymin, periodicity, c1_real, c1_imag, imode, gq_y_points, c1_f)
            call idft_mode_points(ymin, periodicity, c2_real, c2_imag, imode, gq_y_points, c2_f)
            call idft_mode_points(ymin, periodicity, c3_real, c3_imag, imode, gq_y_points, c3_f)
            call idft_mode_points(ymin, periodicity, c4_real, c4_imag, imode, gq_y_points, c4_f)


            !
            ! Compute velocity components
            !
            u_f = rhou_f/rho_f
            v_f = rhov_f/rho_f
            w_f = rhow_f/rho_f


            !
            ! Compute interior speed of sound
            !
            c_f = sqrt(gam_f * p_f / rho_f )


            !
            ! Compute velocity magnitude from interior state
            !
            vmag2_f = (u_f*u_f) + (v_f*v_f) + (w_f*w_f)
            vmag_f = sqrt(vmag2_f)


            !
            ! Compute face total enthalpy and R+ characteristic
            !
            Ht_f    = (p_f / rho_f) * (gam_f/(gam_f - ONE)) + HALF*(vmag2_f)
            Rplus_f = -vmag_f - TWO*c_f/(gam_f - ONE)


            !
            ! solve quadratic equation for c_b
            !
            a = ONE + TWO/(gam_f - ONE)
            b = TWO * Rplus_f
            c = ((gam_f - ONE)/TWO) * (Rplus_f**TWO - TWO*Ht_f)
            

            cb_plus  = (-b  +  sqrt(b**TWO - FOUR*a*c))/(TWO*a)
            cb_minus = (-b  -  sqrt(b**TWO - FOUR*a*c))/(TWO*a)

            c_b = max(cb_plus,cb_minus)



            !
            ! Get boundary condition velocity from extrapolated characteristic and computed c_b
            !
            vmag_b = -TWO*c_b/(gam_m - ONE) - Rplus_f


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
            ! Evaluate influence of characteristic harmonics.
            !
            ! Dont use mean component here so only use 2,nmodes.
            !
            do imode = 2,nmodes

                !
                ! Evaluate modes of characteristics at gq points. c1,c2,c3=0 ( 1D, unsteady, nonreflecting).
                !
                dc4 = idft_mode_points(c4_modes,ymin,periodicity,gq_y_points)

                !
                ! Compute update in primitive variables resulting from characteristics
                !
                rho_b = rho_b + dc4*(ONE/(TWO*(cbar**TWO)))
                u_b   = u_b   - dc4*(ONE/(TWO*rhobar*cbar))
                v_b   = v_b
                w_b   = w_b
                p_b   = p_b   + dc4*HALF
                

            end do ! imode





































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

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irho,integrand)

            !=================================================
            ! x-momentum flux
            !=================================================
            flux_x = (rho_b * u_b * u_b) + p_b
            flux_y = (rho_b * u_b * v_b)
            flux_z = (rho_b * u_b * w_b)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irhou,integrand)

            !=================================================
            ! y-momentum flux
            !=================================================
            flux_x = (rho_b * v_b * u_b)
            flux_y = (rho_b * v_b * v_b) + p_b
            flux_z = (rho_b * v_b * w_b)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irhov,integrand)

            !=================================================
            ! z-momentum flux
            !=================================================
            flux_x = (rho_b * w_b * u_b)
            flux_y = (rho_b * w_b * v_b)
            flux_z = (rho_b * w_b * w_b) + p_b

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irhow,integrand)


            !=================================================
            ! Energy flux
            !=================================================
            flux_x = (rho_b * u_b * H_b)
            flux_y = (rho_b * v_b * H_b)
            flux_z = (rho_b * w_b * H_b)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irhoE,integrand)


        end associate

    end subroutine compute
    !*********************************************************************************************






end module bc_euler_giles_inlet
