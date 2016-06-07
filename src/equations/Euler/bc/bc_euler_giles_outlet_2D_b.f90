module bc_euler_giles_outlet_2D_b
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
    type, public, extends(bc_t) :: euler_giles_outlet_2D_b_t

        type(point_t),  allocatable     :: dft_points(:)

    contains

        procedure   :: add_options              !< Add boundary condition options
        procedure   :: init_spec                !< Specialized bc initialization.
        procedure   :: init_boundary_coupling   !< Implement specialized coupling information between elements.
        procedure   :: compute                  !< boundary condition function implementation

    end type euler_giles_outlet_2D_b_t
    !****************************************************************************************






contains




    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/20/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_options(self)    
        class(euler_giles_outlet_2D_b_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('euler_giles_outlet_2D_b')


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
    subroutine init_spec(self,mesh)
        class(euler_giles_outlet_2D_b_t), intent(inout)   :: self
        type(mesh_t),                intent(inout)   :: mesh

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



        print*, 'dft points'
        do ipt = 1,size(self%dft_points)
            print*, self%dft_points(ipt)%c1_, self%dft_points(ipt)%c2_, self%dft_points(ipt)%c3_
        end do


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
        class(euler_giles_outlet_2D_b_t),   intent(inout)   :: self
        type(mesh_t),                       intent(in)      :: mesh

        integer(ik) :: ielem_bc, ielem_coupled, ielem, ielem_test, var, mode, i

        logical :: same_span = .false.


        print*, 'init_boundary_coupling - begin'

        !
        ! Loop through elements. For the current 2D giles, every element on the boundary
        ! is coupled with every other element on the boundary through the Fourier transform.
        !
        do ielem_bc = 1,size(self%elems)

            print*, 'ielem_bc', ielem_bc 

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

                
                !var  = 1    ! x-span
                !var  = 2    ! y-span
                var  = 3    ! z-span
                mode = 1
                
                same_span = ( abs(mesh%elems(ielem)%coords%getterm(var,mode) - mesh%elems(ielem_test)%coords%getterm(var,mode)) < 0.00001_rk )
                                

                ! Set up boundary-global coupling
                if ( same_span ) then

                    !
                    ! Add element index to the coupling for the current element.
                    !
                    call self%coupled_elems(ielem_bc)%push_back(ielem_test)

                    print*, 'coupled with: ', ielem_test

                end if ! same_span






            end do ! ielem_coupled



        end do  !ielem_bc




     


        print*, 'init_boundary_coupling - end'

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
        class(euler_giles_outlet_2D_b_t),    intent(inout)   :: self
        type(mesh_t),                   intent(in)      :: mesh(:)
        type(solverdata_t),             intent(inout)   :: sdata
        class(properties_t),            intent(inout)   :: prop
        type(face_info_t),              intent(in)      :: face
        type(function_info_t),          intent(in)      :: flux


        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE

        integer(ik)     :: ipt, imode, nmodes

        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(face%idomain)%faces(face%ielement,face%iface)%gq%face%nnodes)   ::  &
                        rho_m,   rhou_m,   rhov_m,   rhow_m,   rhoE_m,          &
                        p_m,     gam_m,    u_m,      v_m,      w_m,             &
                        rho_bar_gq, rhou_bar_gq, rhov_bar_gq, rhow_bar_gq, rhoE_bar_gq,        & 
                        u_bar_gq,   v_bar_gq,    w_bar_gq,    c_bar_gq,    p_bar_gq,  gam_bar_gq, &
                        flux_x,  flux_y,   flux_z,   integrand, tmp,            &
                        c4_mean,                                                &
                        drho_gq, du_gq, dv_gq, dw_gq, dp_gq,                    &
                        drho_bc, du_bc, dv_bc, dp_bc,                           &
                        drho_mode, du_mode, dv_mode, dw_mode, dp_mode,          &
                        drho_mean, du_mean, dv_mean, dw_mean, dp_mean,          &
                        rho_bc, u_bc, v_bc, w_bc, p_bc, rhoE_bc, H_bc,          &
                        B, A2_real, A2_imag, A3_real, A3_imag,                  &
                        c1_gq,    c2_gq,     c3_gq,     c4_gq

        type(AD_D), dimension(size(self%dft_points))    ::                      &
                        rho_b,  rhou_b, rhov_b, rhow_b, rhoE_b,                 &
                        drho,   du,     dv,     dw,     dp,                     &
                        u_b,    v_b,    w_b,    c_b,    p_b,    gam_b,          &
                        rho_bar_b, rhou_bar_b, rhov_bar_b, rhow_bar_b, rhoE_bar_b,        & 
                        u_bar_b,   v_bar_b,    w_bar_b,    c_bar_b,    p_bar_b,  gam_bar_b, &
                        c1_b,   c2_b,   c3_b,   c4_b
                    

        

        ! DFT values
        type(AD_D), allocatable, dimension(:)   ::                              &
                        rho_real, rhou_real, rhov_real, rhow_real, rhoE_real,   &
                        rho_imag, rhou_imag, rhov_imag, rhow_imag, rhoE_imag,   &
                        c1_real,  c2_real,   c3_real,   c4_real,                &
                        c1_imag,  c2_imag,   c3_imag,   c4_imag,                &
                        ptest_real, ptest_imag


        real(rk)        :: periodicity
        real(rk)        :: zero_time
        type(point_t)   :: zero_point

        ! User-specified average pressure
        real(rk),   dimension(mesh(face%idomain)%faces(face%ielement,face%iface)%gq%face%nnodes)   :: p_set, gq_points_y

        real(rk)                :: zavg


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
                ! Get boundary periodicity from bc options.
                !
                zero_time  = ZERO
                call zero_point%set(ZERO,ZERO,ZERO)
                periodicity = self%bcproperties%compute("periodicity", zero_time, zero_point)



                !
                ! Set span location for interpolating on the boundary.
                !
                ! TODO: Hardcoded zavg
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
                ! Compute primitive variables across boundary
                !
                u_b = rhou_b / rho_b
                v_b = rhov_b / rho_b
                w_b = rhow_b / rho_b
                call prop%fluid%compute_pressure(rho_b,rhou_b,rhov_b,rhow_b,rhoE_b,p_b)



                !
                ! DFT conservative variables
                !
                call dft(rho_b,  rho_real,  rho_imag )
                call dft(rhou_b, rhou_real, rhou_imag)
                call dft(rhov_b, rhov_real, rhov_imag)
                call dft(rhow_b, rhow_real, rhow_imag)
                call dft(rhoE_b, rhoE_real, rhoE_imag)



                !
                ! Get boundary mean components at gq nodes
                !
                rho_bar_gq  = rho_real(1)
                rhou_bar_gq = rhou_real(1)
                rhov_bar_gq = rhov_real(1)
                rhow_bar_gq = rhow_real(1)
                rhoE_bar_gq = rhoE_real(1)

                u_bar_gq = rhou_bar_gq/rho_bar_gq
                v_bar_gq = rhov_bar_gq/rho_bar_gq
                w_bar_gq = rhow_bar_gq/rho_bar_gq
                call prop%fluid%compute_pressure(rho_bar_gq,rhou_bar_gq,rhov_bar_gq,rhow_bar_gq,rhoE_bar_gq,p_bar_gq)
                call prop%fluid%compute_gamma(rho_bar_gq,rhou_bar_gq,rhov_bar_gq,rhow_bar_gq,rhoE_bar_gq,gam_bar_gq)
                c_bar_gq = sqrt(gam_bar_gq * p_bar_gq / rho_bar_gq )



                !
                ! Get boundary mean components at dft nodes
                !
                rho_bar_b  = rho_real(1)
                rhou_bar_b = rhou_real(1)
                rhov_bar_b = rhov_real(1)
                rhow_bar_b = rhow_real(1)
                rhoE_bar_b = rhoE_real(1)

                u_bar_b = rhou_bar_b/rho_bar_b
                v_bar_b = rhov_bar_b/rho_bar_b
                w_bar_b = rhow_bar_b/rho_bar_b
                call prop%fluid%compute_pressure(rho_bar_b,rhou_bar_b,rhov_bar_b,rhow_bar_b,rhoE_bar_b,p_bar_b)
                call prop%fluid%compute_gamma(rho_bar_b,rhou_bar_b,rhov_bar_b,rhow_bar_b,rhoE_bar_b,gam_bar_b)
                c_bar_b = sqrt(gam_bar_b * p_bar_b / rho_bar_b )




                !
                ! Compute c4 characteristic due to update in mean pressure
                !
                c4_mean = -TWO*(p_bar_gq - p_set)

                drho_mean =  (ONE/(TWO*c_bar_gq**TWO))       * c4_mean
                du_mean   = -(ONE/(TWO*rho_bar_gq*c_bar_gq)) * c4_mean
                dp_mean   =  HALF                            * c4_mean




                !
                ! Compute total perturbation in primitive variables across the boundary
                !
                drho = rho_b - rho_bar_b
                du   = u_b   - u_bar_b
                dv   = v_b   - v_bar_b
                dw   = w_b   - w_bar_b
                dp   = p_b   - p_bar_b


                !
                ! Compute characteristic variables across the boundary, to be DFT'd.
                !
                c1_b =   (-c_bar_b**TWO)*drho     +         ZERO             +         ZERO             +         dp
                c2_b =            ZERO            +         ZERO             +  (rho_bar_b*c_bar_b)*dv  +        ZERO
                c3_b =            ZERO            +  (rho_bar_b*c_bar_b)*du  +         ZERO             +         dp
                c4_b =            ZERO            -  (rho_bar_b*c_bar_b)*du  +         ZERO             +         dp


                
                !
                ! Compute modes of characteristic variables.
                !
                call dft(c1_b,  c1_real,  c1_imag)
                call dft(c2_b,  c2_real,  c2_imag)
                call dft(c3_b,  c3_real,  c3_imag)
                call dft(c4_b,  c4_real,  c4_imag)



                !
                ! Get y-component of cartesian coordinates for quadrature nodes to evaluate the DFT modes.
                !
                gq_points_y = coords(:)%c2_




                !
                ! initialize derivative arrays
                !
                drho_mode = c_bar_gq*ZERO
                du_mode   = c_bar_gq*ZERO
                dv_mode   = c_bar_gq*ZERO
                dp_mode   = c_bar_gq*ZERO
                c1_gq     = c_bar_gq*ZERO
                c2_gq     = c_bar_gq*ZERO
                c3_gq     = c_bar_gq*ZERO
                c4_gq     = c_bar_gq*ZERO



                !
                ! Loop through Fourier modes and apply boundary conditions.
                !
                nmodes = size(rho_real)
                do imode = 2,(nmodes-1)/2
                !do imode = 1,(nmodes-1)/2

                    !
                    ! 2D steady Giles
                    !
                    B = sqrt( c_bar_gq**TWO - u_bar_gq**TWO - v_bar_gq**TWO )
                    A2_real = - TWO * u_bar_gq * v_bar_gq / ( v_bar_gq**TWO + B**TWO )
                    A2_imag = - TWO * u_bar_gq * B / ( v_bar_gq**TWO + B**TWO )

                    A3_real = ( v_bar_gq**TWO - B**TWO ) / ( v_bar_gq**TWO + B**TWO )
                    A3_imag = ( TWO * v_bar_gq * B ) / ( v_bar_gq**TWO + B**TWO )




                    !
                    ! Apply boundary condition to 4th characteristic of current mode.
                    !
                    c4_real(imode) = A2_real(1)*c2_real(imode)  -  A2_imag(1)*c2_imag(imode)  +  A3_real(1)*c3_real(imode)  -  A3_imag(1)*c3_imag(imode)
                    c4_imag(imode) = A2_real(1)*c2_imag(imode)  +  A2_imag(1)*c2_real(imode)  +  A3_real(1)*c3_imag(imode)  +  A3_imag(1)*c3_real(imode)




                    !
                    ! Get primitive variable perturbations from DFT mode at GQ nodes
                    !
                    !   TODO: HARDCODED -ONE as ymin here.
                    !
                    !call idft_mode_points(-ONE,periodicity, c1_real(imode), c1_imag(imode), imode, gq_points_y, c1_gq)
                    !call idft_mode_points(-ONE,periodicity, c2_real(imode), c2_imag(imode), imode, gq_points_y, c2_gq)
                    !call idft_mode_points(-ONE,periodicity, c3_real(imode), c3_imag(imode), imode, gq_points_y, c3_gq)
                    !call idft_mode_points(-ONE,periodicity, c4_real(imode), c4_imag(imode), imode, gq_points_y, c4_gq)

                    !call idft_mode_points(-0.08205_rk,periodicity, c1_real(imode), c1_imag(imode), imode, gq_points_y, c1_gq)
                    !call idft_mode_points(-0.08205_rk,periodicity, c2_real(imode), c2_imag(imode), imode, gq_points_y, c2_gq)
                    !call idft_mode_points(-0.08205_rk,periodicity, c3_real(imode), c3_imag(imode), imode, gq_points_y, c3_gq)
                    !call idft_mode_points(-0.08205_rk,periodicity, c4_real(imode), c4_imag(imode), imode, gq_points_y, c4_gq)

                    call idft_mode_points(-0._rk,periodicity, c1_real(imode), c1_imag(imode), imode, gq_points_y, c1_gq)
                    call idft_mode_points(-0._rk,periodicity, c2_real(imode), c2_imag(imode), imode, gq_points_y, c2_gq)
                    call idft_mode_points(-0._rk,periodicity, c3_real(imode), c3_imag(imode), imode, gq_points_y, c3_gq)
                    call idft_mode_points(-0._rk,periodicity, c4_real(imode), c4_imag(imode), imode, gq_points_y, c4_gq)

                    !
                    ! Add contribution of current mode characteristics to primitive variable perturbation
                    !
                    drho_mode = drho_mode +     (ONE/(TWO*c_bar_gq**TWO))*c4_gq
                    du_mode   = du_mode   - (ONE/(TWO*rho_bar_gq*c_bar_gq))*c4_gq
                    dv_mode   = dv_mode   +                ZERO
                    dp_mode   = dp_mode   +             HALF*c4_gq



                end do ! imode


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
                drho_gq = rho_m - rho_bar_gq
                du_gq   = u_m   - u_bar_gq
                dv_gq   = v_m   - v_bar_gq
                dw_gq   = w_m   - w_bar_gq
                dp_gq   = p_m   - p_bar_gq


                !
                ! Compute c1,c2,c3 characteristics
                !
                !
                ! Compute characteristic variables across the boundary
                !
                c1_gq =   (-c_bar_gq**TWO)*drho_gq  +              ZERO             +              ZERO             +         dp_gq
                c2_gq =            ZERO             +              ZERO             +  (rho_bar_gq*c_bar_gq)*dv_gq  +        ZERO
                c3_gq =            ZERO             +  (rho_bar_gq*c_bar_gq)*du_gq  +              ZERO             +         dp_gq


                !
                ! Get contribution to primitive variables from c1,c2,c3
                !
                drho_mode = drho_mode + (-ONE/(c_bar_gq**TWO))* c1_gq +              ZERO                +    (ONE/(TWO*c_bar_gq**TWO))    * c3_gq 
                du_mode   = du_mode   +         ZERO                  +              ZERO                + (ONE/(TWO*rho_bar_gq*c_bar_gq)) * c3_gq
                dv_mode   = dv_mode   +         ZERO                  + ONE/(rho_bar_gq*c_bar_gq)* c2_gq +                ZERO                
                dp_mode   = dp_mode   +         ZERO                  +              ZERO                +                HALF             * c3_gq





                !
                ! Compute total primitive variables from mean, contribution from modal update, contribution from mean update.
                !
                rho_bc = rho_bar_gq + drho_mean + drho_mode
                u_bc   = u_bar_gq   + du_mean   + du_mode
                v_bc   = v_bar_gq               + dv_mode
                w_bc   = w_bar_gq
                p_bc   = p_bar_gq   + dp_mean   + dp_mode



                !
                ! Compute boundary condition energy and enthalpy
                !
                rhoE_bc = p_bc/(gam_bar_gq - ONE) + (rho_bc/TWO)*(u_bc*u_bc + v_bc*v_bc + w_bc*w_bc)
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



















end module bc_euler_giles_outlet_2D_b
