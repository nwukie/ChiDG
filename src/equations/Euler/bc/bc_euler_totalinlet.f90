module bc_euler_totalinlet
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ONE, TWO, HALF, ZERO, ME

    use type_bc,            only: bc_t
    use type_solverdata,    only: solverdata_t
    use type_mesh,          only: mesh_t
    use type_properties,    only: properties_t
    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t


    use mod_integrate,      only: integrate_boundary_scalar_flux
    use mod_interpolate,    only: interpolate
    use mod_interpolation,  only: interpolate
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
    !*******************************************************************************************




contains




    !>  Add options for total pressure/temperature boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
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

        call self%bcproperties%add('normal_direction','Required')
        call self%bcproperties%add('nx',              'Required')
        call self%bcproperties%add('ny',              'Required')
        call self%bcproperties%add('nz',              'Required')
        call self%bcproperties%add('nr',              'Required')
        call self%bcproperties%add('nt',              'Required')

        !
        ! Set default angle
        !
        call self%set_fcn_option('nx', 'val', 1._rk)
        call self%set_fcn_option('ny', 'val', 0._rk)
        call self%set_fcn_option('nz', 'val', 0._rk)


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
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,face,fcn)
        class(euler_totalinlet_t),      intent(inout)   :: self
        type(mesh_t),                   intent(in)      :: mesh(:)
        type(solverdata_t),             intent(inout)   :: sdata
        class(properties_t),            intent(inout)   :: prop
        type(face_info_t),              intent(in)      :: face
        type(function_info_t),          intent(in)      :: fcn


        ! Equation indices
        integer(ik)     :: irho, irhou, irhov, irhow, irhoE


        ! Storage at quadrature nodes
        !type(AD_D), dimension(mesh(face%idomain_l)%faces(face%ielement_l,face%iface)%gq%face%nnodes)   ::  &
        type(AD_D), allocatable, dimension(:)   ::  &
                        rho_m,  rhou_m, rhov_m, rhow_m, rhoE_m, p_m,        &
                        flux_x, flux_y, flux_z, integrand,                  &
                        u_m,    v_m,    w_m,                                &
                        u_bc,   v_bc,   w_bc,                               &
                        T_bc,   p_bc,   rho_bc, rhoE_bc,                    &
                        vmag2_m, vmag, H_bc

        real(rk), dimension(mesh(face%idomain_l)%faces(face%ielement_l,face%iface)%gq%face%nnodes)     ::  &
                        TT, PT, nx, ny, nz, nr, nt, theta, normal_direction, x, y, r


        real(rk),   dimension(:), allocatable   :: nt_list, nz_list, r_list

        real(rk)        :: gam_m, cp_m, M
        !real(rk)        :: norm_bc(3)
        integer(ik)     :: iface_p, ineighbor, idonor, igq
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


        idom  = face%idomain_l
        ielem = face%ielement_l
        iface = face%iface



        associate (norms => mesh(idom)%faces(ielem,iface)%norm, unorms => mesh(idom)%faces(ielem,iface)%unorm, faces => mesh(idom)%faces, &
                coords => mesh(idom)%faces(ielem,iface)%quad_pts,    q => sdata%q,      time => sdata%t )


            !
            ! Get boundary condition Total Temperature, Total Pressure, and normal vector
            !
            PT = self%bcproperties%compute("TotalPressure",     time, coords)
            TT = self%bcproperties%compute("TotalTemperature",  time, coords)




            normal_direction = self%bcproperties%compute("normal_direction",    time, coords)

!            if ( normal_direction(1) == ONE ) then
                nx = self%bcproperties%compute("nx",                                time, coords)
                ny = self%bcproperties%compute("ny",                                time, coords)
                nz = self%bcproperties%compute("nz",                                time, coords)

!            else if ( normal_direction(1) == TWO ) then
!                nr = self%bcproperties%compute("nr",                                time, coords)
!                nt = self%bcproperties%compute("nt",                                time, coords)
!                nz = self%bcproperties%compute("nz",                                time, coords)
!
!
!                !
!                ! Compute theta
!                !
!                x = mesh(idom)%elems(ielem)%quad_pts(:)%c1_
!                y = mesh(idom)%elems(ielem)%quad_pts(:)%c2_
!
!                theta = atan2(y,x)
!
!                !
!                ! Convert cylindrical vectors to cartesian
!                ! 
!                nx = cos(theta)*nr - sin(theta)*nt
!                ny = sin(theta)*nr + cos(theta)*nt
!                nz = nz
!
!            else
!                call write_line("Warning: bc_euler_totalinlet not reading inlet direction")
!
!            end if


            



            !
            ! Interpolate interior solution to quadrature nodes
            !
!            call interpolate_face(mesh,face,fcn,q,irho, rho_m,  'value', ME)
!            call interpolate_face(mesh,face,fcn,q,irhou,rhou_m, 'value', ME)
!            call interpolate_face(mesh,face,fcn,q,irhov,rhov_m, 'value', ME)
!            call interpolate_face(mesh,face,fcn,q,irhow,rhow_m, 'value', ME)
!            call interpolate_face(mesh,face,fcn,q,irhoE,rhoE_m, 'value', ME)

            rho_m  = interpolate(mesh,sdata,face,fcn,irho,  'value', ME)
            rhou_m = interpolate(mesh,sdata,face,fcn,irhou, 'value', ME)
            rhov_m = interpolate(mesh,sdata,face,fcn,irhov, 'value', ME)
            rhow_m = interpolate(mesh,sdata,face,fcn,irhow, 'value', ME)
            rhoE_m = interpolate(mesh,sdata,face,fcn,irhoE, 'value', ME)


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
            u_bc = vmag*nx
            v_bc = vmag*ny
            w_bc = vmag*nz








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

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irho,integrand)

            !=================================================
            ! x-momentum flux
            !=================================================
            flux_x = (rho_bc * u_bc * u_bc) + p_bc
            flux_y = (rho_bc * u_bc * v_bc)
            flux_z = (rho_bc * u_bc * w_bc)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irhou,integrand)

            !=================================================
            ! y-momentum flux
            !=================================================
            flux_x = (rho_bc * v_bc * u_bc)
            flux_y = (rho_bc * v_bc * v_bc) + p_bc
            flux_z = (rho_bc * v_bc * w_bc)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irhov,integrand)

            !=================================================
            ! z-momentum flux
            !=================================================
            flux_x = (rho_bc * w_bc * u_bc)
            flux_y = (rho_bc * w_bc * v_bc)
            flux_z = (rho_bc * w_bc * w_bc) + p_bc

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irhow,integrand)


            !=================================================
            ! Energy flux
            !=================================================
            flux_x = (rho_bc * u_bc * H_bc)
            flux_y = (rho_bc * v_bc * H_bc)
            flux_z = (rho_bc * w_bc * H_bc)

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irhoE,integrand)


        end associate

    end subroutine compute
    !***********************************************************************************************






end module bc_euler_totalinlet
