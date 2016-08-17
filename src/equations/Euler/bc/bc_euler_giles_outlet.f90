module bc_euler_giles_outlet
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO, ONE, TWO, HALF, ME
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
    type, public, extends(bc_t) :: euler_giles_outlet_t

        type(point_t),  allocatable     :: dft_points(:)
!        real(rk),       allocatable     :: dft_span_z(:)

    contains

        procedure   :: add_options              !< Add boundary condition options
        procedure   :: init_spec                !< Specialized bc initialization.
        procedure   :: init_boundary_coupling   !< Implement specialized coupling information between elements.
        procedure   :: compute                  !< boundary condition function implementation

    end type euler_giles_outlet_t
    !****************************************************************************************




contains


    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/20/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_options(self)    
        class(euler_giles_outlet_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('euler_giles_outlet')


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
        class(euler_giles_outlet_t), intent(inout)   :: self
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


!        !
!        ! Compute span locations along the boundary. Assuming span=z.
!        !
!        self%dft_span_z = compute_span_locations(mesh,self%elems,iface,periodicity)


        !
        ! Compute dft points
        !
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
        class(euler_giles_outlet_t), intent(inout)   :: self
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
    subroutine compute(self,mesh,sdata,prop,face,fcn)
        class(euler_giles_outlet_t),    intent(inout)   :: self
        type(mesh_t),                   intent(in)      :: mesh(:)
        type(solverdata_t),             intent(inout)   :: sdata
        class(properties_t),            intent(inout)   :: prop
        type(face_info_t),              intent(in)      :: face
        type(function_info_t),          intent(in)      :: fcn


        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE

        integer(ik)     :: ipt

        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(face%idomain)%faces(face%ielement,face%iface)%gq%face%nnodes)   ::  &
                        rho_m,  rhou_m, rhov_m, rhow_m, rhoE_m,             &
                        p_m,    gam_m,  u_m,    v_m,    w_m,                &
                        p_bc,   H_bc,   rhoE_bc,                            &
                        rho_bar, u_bar, c_bar,                              &
                        flux_x, flux_y, flux_z, integrand, tmp, c3

        type(AD_D), dimension(size(self%dft_points))    :: &
                        rho_b,  rhou_b, rhov_b, rhow_b, rhoE_b,             &
                        u_b,    c_b,    p_b,    gam_b
                    


        ! DFT values
        type(AD_D), allocatable, dimension(:)   :: rho_real, rho_imag, u_real, u_imag, c_real, c_imag

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



                call interpolate_boundary(mesh,face,q,irho , self%dft_points, rho_b )
                call interpolate_boundary(mesh,face,q,irhou, self%dft_points, rhou_b)
                call interpolate_boundary(mesh,face,q,irhov, self%dft_points, rhov_b)
                call interpolate_boundary(mesh,face,q,irhow, self%dft_points, rhow_b)
                call interpolate_boundary(mesh,face,q,irhoE, self%dft_points, rhoE_b)
                



                !
                ! Compute u-velocity and pressure across boundary. To be DFT'd.
                !
                call prop%fluid%compute_pressure(rho_b,rhou_b,rhov_b,rhow_b,rhoE_b,p_b)
                u_b = rhou_b / rho_b



                !
                ! Compute speed of sound across boundary
                !
                call prop%fluid%compute_gamma(rho_b,rhou_b,rhov_b,rhow_b,rhoE_b,gam_b)
                c_b = sqrt(gam_b * p_b / rho_b )



                !
                ! DFT rho, u, c to get average values across face. average p is set by user: p_set.
                !
                call dft(rho_b,  rho_real,  rho_imag)
                call dft(u_b,    u_real,    u_imag  )
                call dft(c_b,    c_real,    c_imag  )




                
                !
                ! Get DC component from DFT.
                !
                rho_bar = rho_real(1)
                u_bar   = u_real(1)
                c_bar   = c_real(1)



                !
                ! Interpolate interior solution to face quadrature nodes
                !
                call interpolate_face(mesh,face,fcn,q,irho, rho_m,  'value', ME)
                call interpolate_face(mesh,face,fcn,q,irhou,rhou_m, 'value', ME)
                call interpolate_face(mesh,face,fcn,q,irhov,rhov_m, 'value', ME)
                call interpolate_face(mesh,face,fcn,q,irhow,rhow_m, 'value', ME)
                call interpolate_face(mesh,face,fcn,q,irhoE,rhoE_m, 'value', ME)



                !
                ! Compute pressure, gamma from interior.
                !
                call prop%fluid%compute_pressure(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,p_m  )
                call prop%fluid%compute_gamma(   rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,gam_m)



                !
                ! Compute velocity components
                !
                u_m = rhou_m/rho_m
                v_m = rhov_m/rho_m
                w_m = rhow_m/rho_m




                !
                ! Compute pressure. IMPOSED BOUNDARY CONDITION.
                !
                p_bc = HALF*( (rho_bar*c_bar)*(u_m - u_bar) + (p_m + p_set) )


                !c3 = (rho_bar * c_bar)*(u_m - u_bar) + (p_m - p_set)
                !p_bc = p_set + HALF*c3




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


                call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irho,integrand)


                !=================================================
                ! x-momentum flux
                !=================================================
                flux_x = (rho_m * u_m * u_m) + p_bc
                flux_y = (rho_m * u_m * v_m)
                flux_z = (rho_m * u_m * w_m)


                integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)


                call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irhou,integrand)


                !=================================================
                ! y-momentum flux
                !=================================================
                flux_x = (rho_m * v_m * u_m)
                flux_y = (rho_m * v_m * v_m) + p_bc
                flux_z = (rho_m * v_m * w_m)

                integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

                call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irhov,integrand)


                !=================================================
                ! z-momentum flux
                !=================================================
                flux_x = (rho_m * w_m * u_m)
                flux_y = (rho_m * w_m * v_m)
                flux_z = (rho_m * w_m * w_m) + p_bc

                integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

                call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irhow,integrand)


                !=================================================
                ! Energy flux
                !=================================================
                flux_x = (rho_m * u_m * H_bc)
                flux_y = (rho_m * v_m * H_bc)
                flux_z = (rho_m * w_m * H_bc)

                integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

                call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irhoE,integrand)


            end associate

        end associate






    end subroutine compute
    !**********************************************************************************************












!    !>
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   2/21/2016
!    !!
!    !!
!    !!
!    !!
!    !-----------------------------------------------------------------------------------------------
!    function compute_span_locations(mesh,elems) result(spans)
!        type(mesh_t),   intent(in)  :: mesh
!        integer(ik),    intent(in)  :: elems(:)
!
!        type(rvector_t) :: span_locations
!        integer(ik)     :: var, term, ielem_bc, ielem, ispan
!        real(rk)        :: zavg
!        logical         :: span_already_exists
!
!        real(rk),   allocatable :: spans(:)
!
!        !
!        ! Compute number of unique span stations.
!        !
!        do ielem_bc = 1,size(elems)
!
!            ielem = elems(ielem_bc)
!            
!            !
!            ! Get average z of current element
!            !
!            var  = 3
!            term = 1
!            zavg = mesh%elems(ielem)%coords%getterm(3,1)
!
!
!            !
!            ! Check if zavg already exists in the recorded span locations.
!            !
!            span_already_exists = .false.
!            do ispan = 1,span_locations%size()
!
!                iz = span_locations%at(ispan)
!
!                span_already_exists = ( abs( iz - zavg ) < 0.01_rk )
!
!                if ( span_already_exists ) then
!                    exit
!                end if
!
!            end do ! ispan
!
!            !
!            ! If current zavg was not detected in the list, store as new span location.
!            !
!            if ( .not. span_already_exists ) then
!                call span_locations%push_back(zavg)
!            end if
!
!        end do ! ielem
!
!
!
!
!        !
!        ! Assemble array to return
!        !
!        nspan = span_locations%size()
!        allocate( spans(nspan), stat=ierr)
!        if (ierr /= 0) call AllocationError
!
!        do ispan = 1,nspan
!            spans(ispan) = span_locations%at(ispan)
!        end do
!
!
!    end function compute_span_locations
!    !************************************************************************************************
!
















end module bc_euler_giles_outlet
