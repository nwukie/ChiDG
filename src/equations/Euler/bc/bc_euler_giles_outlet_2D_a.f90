module bc_euler_giles_outlet_2D_a
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO, ONE, TWO, HALF, LOCAL
    use type_bc,            only: bc_t
    use type_solverdata,    only: solverdata_t
    use type_point,         only: point_t
    use type_mesh,          only: mesh_t
    use type_properties,    only: properties_t
    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t

    use mod_integrate,      only: integrate_boundary_scalar_flux
    use mod_interpolate,    only: interpolate_face, interpolate_boundary
    use mod_dft,            only: dft, idft_mode_points, compute_dft_points
    use DNAD_D
    
    use EULER_properties,   only: EULER_properties_t
    implicit none


    !> Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   2/20/2016
    !!
    !----------------------------------------------------------------------------------------
    type, public, extends(bc_t) :: euler_giles_outlet_2D_a_t

        type(point_t),  allocatable     :: dft_points(:)
!        real(rk),       allocatable     :: dft_span_z(:)

    contains

        procedure   :: add_options              !< Add boundary condition options
        procedure   :: init_spec                !< Specialized bc initialization.
        procedure   :: init_boundary_coupling   !< Implement specialized coupling information between elements.
        procedure   :: compute                  !< boundary condition function implementation

    end type euler_giles_outlet_2D_a_t
    !****************************************************************************************




contains


    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/20/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_options(self)    
        class(euler_giles_outlet_2D_a_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('euler_giles_outlet_2D_a')


        !
        ! Add functions
        !
        call self%bcproperties%add('StaticPressure','Required')         ! add StaticPressure
        call self%bcproperties%add('periodicity',   'Required')


        !
        ! Add parameters
        !


    end subroutine add_options
    !******************************************************************************************







    !>  Specialized initialization for the boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/16/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------
    !subroutine init_spec(self,mesh,iface)
    subroutine init_spec(self,mesh)
        class(euler_giles_outlet_2D_a_t), intent(inout)   :: self
        type(mesh_t),                intent(inout)   :: mesh
        !integer(ik),                 intent(in)      :: iface

        real(rk)        :: periodicity
        real(rk)        :: zero_time
        type(point_t)   :: zero_point

        integer(ik) :: ipt

        !
        ! Get boundary periodicity from bc options.
        !
        zero_time  = ZERO
        call zero_point%set(ZERO,ZERO,ZERO)
        periodicity = self%bcproperties%compute("periodicity", zero_time, zero_point)



        !
        ! Compute dft points
        !
        !self%dft_points = compute_dft_points(mesh,self%elems,iface,periodicity)
        self%dft_points = compute_dft_points(mesh,self%elems,self%faces,periodicity)


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
    subroutine init_boundary_coupling(self,mesh)
        class(euler_giles_outlet_2D_a_t), intent(inout)   :: self
        type(mesh_t),                intent(in)      :: mesh

        integer(ik) :: ielem_bc, ielem_coupled, ielem, ielem_test, var, mode, i

        logical :: same_span = .false.

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
                ! Get block-element index of potentially coupled element
                !
                ielem_test = self%elems(ielem_coupled)

                
                var  = 3
                mode = 1
                
                same_span = ( abs(mesh%elems(ielem)%coords%getterm(3,mode) - mesh%elems(ielem_test)%coords%getterm(3,mode)) < 0.00001_rk )
                                

                if ( same_span ) then

                    !
                    ! Add element index to the coupling for the current element.
                    !
                    call self%coupled_elems(ielem_bc)%push_back(ielem_test)

                end if ! same_span




            end do ! ielem_coupled


        end do  !ielem_bc



    end subroutine init_boundary_coupling
    !*********************************************************************************************















    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!  @param[in]      face    face_info_t containing indices on location and misc information
    !!  @param[in]      flux    function_into_t containing info on the function being computed
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,face,flux)
        class(euler_giles_outlet_2D_a_t),    intent(inout)   :: self
        type(mesh_t),                   intent(in)      :: mesh(:)
        type(solverdata_t),             intent(inout)   :: sdata
        class(properties_t),            intent(inout)   :: prop
        type(face_info_t),              intent(in)      :: face
        type(function_info_t),          intent(in)      :: flux


        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE

        integer(ik)     :: ipt

        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(face%idomain)%faces(face%ielement,face%iface)%gq%face%nnodes)   ::  &
                        rho_m,   rhou_m,   rhov_m,   rhow_m,   rhoE_m,          &
                        p_m,     gam_m,    u_m,      v_m,      w_m,             &
                        rho_bar, rhou_bar, rhov_bar, rhow_bar, rhoE_bar,        & 
                        u_bar,   v_bar,    w_bar,    c_bar,    p_bar,  gam_bar, &
                        flux_x,  flux_y,   flux_z,   integrand, tmp,            &
                        c1, c2, c3, c4, c4_incoming, c4_mean,                   &
                        drho, du, dv, dw, dp,                                   &
                        drho_bc, du_bc, dv_bc, dp_bc,                           &
                        rho_bc, u_bc, v_bc, w_bc, p_bc, rhoE_bc, H_bc,          &
                        A2, A4, B, beta

        type(AD_D), dimension(size(self%dft_points))    ::                      &
                        rho_b,  rhou_b, rhov_b, rhow_b, rhoE_b,                 &
                        u_b,    c_b,    p_b,    gam_b
                    


        ! DFT values
        type(AD_D), allocatable, dimension(:)   ::                              &
                        rho_real, rhou_real, rhov_real, rhow_real, rhoE_real,   &
                        rho_imag, rhou_imag, rhov_imag, rhow_imag, rhoE_imag

        ! User-specified average pressure
        real(rk),   dimension(mesh(face%idomain)%faces(face%ielement,face%iface)%gq%face%nnodes)   :: p_set

        real(rk)    :: zavg


        associate ( idom => face%idomain, ielem => face%ielement, iface => face%iface )

            associate ( norms  => mesh(idom)%faces(ielem,iface)%norm,        unorms => mesh(idom)%faces(ielem,iface)%unorm, &
                        coords => mesh(idom)%faces(ielem,iface)%quad_pts,    q => sdata%q,      time => sdata%t )


                
                !
                ! Get equation indices
                !
                irho  = prop%get_eqn_index("rho")
                irhou = prop%get_eqn_index("rhou")
                irhov = prop%get_eqn_index("rhov")
                irhow = prop%get_eqn_index("rhow")
                irhoE = prop%get_eqn_index("rhoE")




                !
                ! Get back pressure from function.
                !
                p_set = self%bcproperties%compute("StaticPressure",time,coords)



                !
                ! Set span location for interpolating on the boundary.
                !
                zavg = mesh(idom)%elems(ielem)%coords%getterm(3,1)
                do ipt = 1,size(self%dft_points)
                    call self%dft_points(ipt)%z(zavg)
                end do



                !
                ! Interpolate solution across the boundary to be DFT'd.
                !
                call interpolate_boundary(mesh,face,q,irho , self%dft_points, rho_b )
                call interpolate_boundary(mesh,face,q,irhou, self%dft_points, rhou_b)
                call interpolate_boundary(mesh,face,q,irhov, self%dft_points, rhov_b)
                call interpolate_boundary(mesh,face,q,irhow, self%dft_points, rhow_b)
                call interpolate_boundary(mesh,face,q,irhoE, self%dft_points, rhoE_b)
                


                !
                ! DFT conservative variables
                !
                call dft(rho_b,  rho_real,  rho_imag )
                call dft(rhou_b, rhou_real, rhou_imag)
                call dft(rhov_b, rhov_real, rhov_imag)
                call dft(rhow_b, rhow_real, rhow_imag)
                call dft(rhoE_b, rhoE_real, rhoE_imag)


                !
                ! Get boundary mean components
                !
                rho_bar  = rho_real(1)
                rhou_bar = rhou_real(1)
                rhov_bar = rhov_real(1)
                rhow_bar = rhow_real(1)
                rhoE_bar = rhoE_real(1)

                u_bar = rhou_bar/rho_bar
                v_bar = rhov_bar/rho_bar
                w_bar = rhow_bar/rho_bar
                call prop%fluid%compute_pressure(rho_bar,rhou_bar,rhov_bar,rhow_bar,rhoE_bar,p_bar)
                call prop%fluid%compute_gamma(rho_bar,rhou_bar,rhov_bar,rhow_bar,rhoE_bar,gam_bar)
                c_bar = sqrt(gam_bar * p_bar / rho_bar )



                !
                ! Compute c4 characteristic due to update in mean pressure
                !
                c4_mean = -TWO*(p_bar - p_set)




                !
                ! Interpolate interior solution to face quadrature nodes
                !
                call interpolate_face(mesh,face,q,irho, rho_m, LOCAL)
                call interpolate_face(mesh,face,q,irhou,rhou_m,LOCAL)
                call interpolate_face(mesh,face,q,irhov,rhov_m,LOCAL)
                call interpolate_face(mesh,face,q,irhow,rhow_m,LOCAL)
                call interpolate_face(mesh,face,q,irhoE,rhoE_m,LOCAL)


                !
                ! Compute interior primitive variables
                !
                u_m = rhou_m / rho_m
                v_m = rhov_m / rho_m
                w_m = rhow_m / rho_m
                call prop%fluid%compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,p_m)


                !
                ! Compute variation in primitive variables from the mean.
                !
                drho = rho_m - rho_bar
                du   = u_m   - u_bar
                dv   = v_m   - v_bar
                dw   = w_m   - w_bar
                dp   = p_m   - p_bar



                !
                ! Compute c1,c2,c3,c4 characteristics due to variation in primitive variables
                !
                c1 = (-c_bar**TWO)*drho         +         ZERO         +         ZERO         +         dp
                c2 =        ZERO                +         ZERO         +  (rho_bar*c_bar)*dv  +        ZERO
                c3 =        ZERO                +  (rho_bar*c_bar)*du  +         ZERO         +         dp
                c4 =        ZERO                -  (rho_bar*c_bar)*du  +         ZERO         +         dp


                !
                ! Impose boundary condition on the characteristics coming from the interior.
                !

                !
                ! 1D unsteady Giles
                !
                !c4 = ZERO


                !
                ! 2D steady Giles
                !
                beta = c_bar**TWO - u_bar**two - v_bar**two
                B = sqrt(beta)

                A2 = TWO*u_bar/(B - v_bar)
                A4 = (B+v_bar)/(B - v_bar)

                c4 = A2*c2  -  A4*c3


                !
                ! Compute the incoming characteristic.
                !
                c4_incoming = c4_mean + c4


                !
                ! Compute primitive variable update from modified characteristics.
                !
                drho_bc = (-ONE/(c_bar**TWO))*c1    +          ZERO            +   (ONE/(TWO*c_bar**TWO))*c3      +   (ONE/(TWO*c_bar**TWO))*c4_incoming
                du_bc   =         ZERO              +          ZERO            +   (ONE/(TWO*rho_bar*c_bar))*c3   -   (ONE/(TWO*rho_bar*c_bar))*c4_incoming
                dv_bc   =         ZERO              +  ONE/(rho_bar*c_bar)*c2  +              ZERO                +                  ZERO
                dp_bc   =         ZERO              +          ZERO            +            HALF*c3               +             HALF*c4_incoming


                !
                ! Compute total primitive variables for boundary condition
                !
                rho_bc = rho_bar + drho_bc
                u_bc   = u_bar   + du_bc
                v_bc   = v_bar   + dv_bc
                w_bc   = w_m
                p_bc   = p_bar   + dp_bc





                !
                ! Compute boundary condition energy and enthalpy
                !
                rhoE_bc = p_bc/(gam_bar - ONE) + (rho_bc/TWO)*(u_bc*u_bc + v_bc*v_bc + w_bc*w_bc)
                H_bc = (rhoE_bc + p_bc)/rho_bc






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

        end associate






    end subroutine compute
    !**********************************************************************************************



















end module bc_euler_giles_outlet_2D_a
