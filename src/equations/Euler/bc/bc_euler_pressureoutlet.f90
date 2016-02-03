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

        procedure   :: set_options  !< Set boundary condition options
        procedure   :: compute      !< boundary condition function implementation

    end type euler_pressureoutlet_t
    !****************************************************************************************




contains


    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine set_options(self)    
        class(euler_pressureoutlet_t),  intent(inout)   :: self


        !
        ! Add functions
        !
        call self%bcfunctions%add('Static Pressure','Required')


        !& DEBUG. This should be removed in general
        call self%bcfunctions%set_fcn(       'Static Pressure', 'constant')
        call self%bcfunctions%set_fcn_option('Static Pressure','val',93000._rk)



        !
        ! Add parameters
        !


    end subroutine set_options
    !******************************************************************************************








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
        class(euler_pressureoutlet_t),  intent(inout)   :: self
        type(mesh_t),                   intent(in)      :: mesh(:)
        type(solverdata_t),             intent(inout)   :: sdata
        class(properties_t),            intent(inout)   :: prop
        type(face_info_t),              intent(in)      :: face
        type(function_info_t),          intent(in)      :: flux


        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE


        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(face%idomain)%faces(face%ielement,face%iface)%gq%face%nnodes)   ::  &
                        rho_m,  rhou_m, rhov_m, rhow_m, rhoE_m,             &
                        flux_x, flux_y, flux_z, integrand,                  &
                        u_m,    v_m,    w_m,                                &
                        H_bc,   rhoE_bc, gam_m

        real(rk),   dimension(mesh(face%idomain)%faces(face%ielement,face%iface)%gq%face%nnodes)   :: p_bc


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
                p_bc = self%bcfunctions%compute("Static Pressure",time,coords)
                !p_bc = 93000._rk
                !p_bc = 107000._rk




                !
                ! Interpolate interior solution to face quadrature nodes
                !
                call interpolate_face(mesh,face,q,irho, rho_m, LOCAL)
                call interpolate_face(mesh,face,q,irhou,rhou_m,LOCAL)
                call interpolate_face(mesh,face,q,irhov,rhov_m,LOCAL)
                call interpolate_face(mesh,face,q,irhow,rhow_m,LOCAL)
                call interpolate_face(mesh,face,q,irhoE,rhoE_m,LOCAL)



                !
                ! Compute gamma
                !
                call prop%fluid%compute_gamma(rho_m,rhou_m,rhov_m,rhow_m,rhoE_m,gam_m)



                !
                ! Compute velocity components
                !
                u_m = rhou_m/rho_m
                v_m = rhov_m/rho_m
                w_m = rhow_m/rho_m



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

        end associate
    end subroutine compute
    !**********************************************************************************************






end module bc_euler_pressureoutlet
