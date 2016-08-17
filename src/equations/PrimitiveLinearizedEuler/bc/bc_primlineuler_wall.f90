module bc_primlineuler_wall
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ONE, TWO, HALF, ZERO, ME
    use type_bc,            only: bc_t
    use type_solverdata,    only: solverdata_t
    use type_mesh,          only: mesh_t
    use type_properties,    only: properties_t
    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t

    use mod_integrate,      only: integrate_boundary_scalar_flux
    use mod_interpolate,    only: interpolate_face
    use DNAD_D
    
    use PRIMLINEULER_properties,   only: PRIMLINEULER_properties_t
    use mod_primitive_linearized_euler
    implicit none


    !>  Extrapolation boundary condition 
    !!      - Extrapolate interior variables to be used for calculating the boundary flux.
    !!  
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !-------------------------------------------------------------------------------------------
    type, public, extends(bc_t) :: primlineuler_wall_t

    contains
    
        procedure   :: add_options
        procedure   :: compute    !> bc implementation

    end type primlineuler_wall_t
    !*******************************************************************************************




contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine add_options(self)
        class(primlineuler_wall_t),    intent(inout)   :: self

        !
        ! Set name
        ! 
        call self%set_name('primlineuler_wall')


    end subroutine add_options
    !*******************************************************************************************









    !> Specialized compute routine for Extrapolation Boundary Condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!  @param[in]      mesh    Mesh data containing elements and faces for the domain
    !!  @param[inout]   sdata   Solver data containing solution vector, rhs, linearization, etc.
    !!  @param[in]      ielem   Index of the element being computed
    !!  @param[in]      iface   Index of the face being computed
    !!  @param[in]      iblk    Index of the linearization block being computed
    !!  @param[inout]   prop    properties_t object containing equations and material_t objects
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,face,fcn)
        class(primlineuler_wall_t),     intent(inout)   :: self
        type(mesh_t),                   intent(in)      :: mesh(:)
        type(solverdata_t),             intent(inout)   :: sdata
        class(properties_t),            intent(inout)   :: prop
        type(face_info_t),              intent(in)      :: face
        type(function_info_t),          intent(in)      :: fcn


        ! Equation indices
        integer(ik)     :: irho_r, iu_r, iv_r, iw_r, ip_r
        integer(ik)     :: irho_i, iu_i, iv_i, iw_i, ip_i


        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(face%idomain_l)%faces(face%ielement_l,face%iface)%gq%face%nnodes)   ::  &
                        rho_r,      u_r,     v_r,     w_r,     p_r,         &
                        rho_i,      u_i,     v_i,     w_i,     p_i,         &
                        c1,         c2,         c3,     c4,     &
                        drho,       du,         dp,                                     &
                        drho_total, du_total,   dp_total,                               &
                        drho_user,  du_user,    dp_user,                                &
                        flux_x, flux_y, flux_z, integrand





        !
        ! Get equation indices
        !
        irho_r = prop%get_eqn_index("rho_r")
        iu_r   = prop%get_eqn_index("u_r")
        iv_r   = prop%get_eqn_index("v_r")
        iw_r   = prop%get_eqn_index("w_r")
        ip_r   = prop%get_eqn_index("p_r")


        irho_i = prop%get_eqn_index("rho_i")
        iu_i   = prop%get_eqn_index("u_i")
        iv_i   = prop%get_eqn_index("v_i")
        iw_i   = prop%get_eqn_index("w_i")
        ip_i   = prop%get_eqn_index("p_i")




        associate ( idom => face%idomain_l, ielem => face%ielement_l, iface => face%iface )

            associate (norms => mesh(idom)%faces(ielem,iface)%norm, unorms => mesh(idom)%faces(ielem,iface)%unorm, faces => mesh(idom)%faces, q => sdata%q)

            !
            ! Interpolate interior solution to quadrature nodes
            !
            call interpolate_face(mesh,face,fcn,q,irho_r, rho_r,  'value', ME)
            call interpolate_face(mesh,face,fcn,q,iu_r,   u_r,    'value', ME)
            call interpolate_face(mesh,face,fcn,q,iv_r,   v_r,    'value', ME)
            call interpolate_face(mesh,face,fcn,q,iw_r,   w_r,    'value', ME)
            call interpolate_face(mesh,face,fcn,q,ip_r,   p_r,    'value', ME)


            call interpolate_face(mesh,face,fcn,q,irho_i, rho_i,  'value', ME)
            call interpolate_face(mesh,face,fcn,q,iu_i,   u_i,    'value', ME)
            call interpolate_face(mesh,face,fcn,q,iv_i,   v_i,    'value', ME)
            call interpolate_face(mesh,face,fcn,q,iw_i,   w_i,    'value', ME)
            call interpolate_face(mesh,face,fcn,q,ip_i,   p_i,    'value', ME)






            !=================================================
            ! Mass flux
            !=================================================
            flux_x = rho_x_rho  * rho_r  + &
!                     rho_x_rhou * rhou_r + &
!                     rho_x_rhov * rhov_r + &
!                     rho_x_rhow * rhow_r + &
                     rho_x_p * p_r
            flux_y = rho_y_rho  * rho_r  + &
!                     rho_y_rhou * rhou_r + &
!                     rho_y_rhov * rhov_r + &
!                     rho_y_rhow * rhow_r + &
                     rho_y_p * p_r
            flux_z = rho_z_rho  * rho_r  + &
!                     rho_z_rhou * rhou_r + &
!                     rho_z_rhov * rhov_r + &
!                     rho_z_rhow * rhow_r + &
                     rho_z_p * p_r

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irho_r,integrand)





            flux_x = rho_x_rho  * rho_i  + &
!                     rho_x_rhou * rhou_i + &
!                     rho_x_rhov * rhov_i + &
!                     rho_x_rhow * rhow_i + &
                     rho_x_p * p_i
            flux_y = rho_y_rho  * rho_i  + &
!                     rho_y_rhou * rhou_i + &
!                     rho_y_rhov * rhov_i + &
!                     rho_y_rhow * rhow_i + &
                     rho_y_p * p_i
            flux_z = rho_z_rho  * rho_i  + &
!                     rho_z_rhou * rhou_i + &
!                     rho_z_rhov * rhov_i + &
!                     rho_z_rhow * rhow_i + &
                     rho_z_p * p_i

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,irho_i,integrand)









            !=================================================
            ! x-momentum flux
            !=================================================
            flux_x = u_x_rho  * rho_r  + &
!                     rhou_x_rhou * rhou_r + &
!                     rhou_x_rhov * rhov_r + &
!                     rhou_x_rhow * rhow_r + &
                     u_x_p * p_r
            flux_y = u_y_rho  * rho_r  + &
!                     rhou_y_rhou * rhou_r + &
!                     rhou_y_rhov * rhov_r + &
!                     rhou_y_rhow * rhow_r + &
                     u_y_p * p_r
            flux_z = u_z_rho  * rho_r  + &
!                     rhou_z_rhou * rhou_r + &
!                     rhou_z_rhov * rhov_r + &
!                     rhou_z_rhow * rhow_r + &
                     u_z_p * p_r

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,iu_r,integrand)





            flux_x = u_x_rho  * rho_i  + &
!                     rhou_x_rhou * rhou_i + &
!                     rhou_x_rhov * rhov_i + &
!                     rhou_x_rhow * rhow_i + &
                     u_x_p * p_i
            flux_y = u_y_rho  * rho_i  + &
!                     rhou_y_rhou * rhou_i + &
!                     rhou_y_rhov * rhov_i + &
!                     rhou_y_rhow * rhow_i + &
                     u_y_p * p_i
            flux_z = u_z_rho  * rho_i  + &
!                     rhou_z_rhou * rhou_i + &
!                     rhou_z_rhov * rhov_i + &
!                     rhou_z_rhow * rhow_i + &
                     u_z_p * p_i

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,iu_i,integrand)







            !=================================================
            ! y-momentum flux
            !=================================================

            flux_x = v_x_rho  * rho_r  + &
!                     rhov_x_rhou * rhou_r + &
!                     rhov_x_rhov * rhov_r + &
!                     rhov_x_rhow * rhow_r + &
                     v_x_p * p_r
            flux_y = v_y_rho  * rho_r  + &
!                     rhov_y_rhou * rhou_r + &
!                     rhov_y_rhov * rhov_r + &
!                     rhov_y_rhow * rhow_r + &
                     v_y_p * p_r
            flux_z = v_z_rho  * rho_r  + &
!                     rhov_z_rhou * rhou_r + &
!                     rhov_z_rhov * rhov_r + &
!                     rhov_z_rhow * rhow_r + &
                     v_z_p * p_r

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,iv_r,integrand)




            flux_x = v_x_rho  * rho_i  + &
!                     rhov_x_rhou * rhou_i + &
!                     rhov_x_rhov * rhov_i + &
!                     rhov_x_rhow * rhow_i + &
                     v_x_p * p_i
            flux_y = v_y_rho  * rho_i  + &
!                     rhov_y_rhou * rhou_i + &
!                     rhov_y_rhov * rhov_i + &
!                     rhov_y_rhow * rhow_i + &
                     v_y_p * p_i
            flux_z = v_z_rho  * rho_i  + &
!                     rhov_z_rhou * rhou_i + &
!                     rhov_z_rhov * rhov_i + &
!                     rhov_z_rhow * rhow_i + &
                     v_z_p * p_i

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,iv_i,integrand)








            !=================================================
            ! z-momentum flux
            !=================================================

            flux_x = w_x_rho  * rho_r  + &
!                     rhow_x_rhou * rhou_r + &
!                     rhow_x_rhov * rhov_r + &
!                     rhow_x_rhow * rhow_r + &
                     w_x_p * p_r
            flux_y = w_y_rho  * rho_r  + &
!                     rhow_y_rhou * rhou_r + &
!                     rhow_y_rhov * rhov_r + &
!                     rhow_y_rhow * rhow_r + &
                     w_y_p * p_r
            flux_z = w_z_rho  * rho_r  + &
!                     rhow_z_rhou * rhou_r + &
!                     rhow_z_rhov * rhov_r + &
!                     rhow_z_rhow * rhow_r + &
                     w_z_p * p_r

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,iw_r,integrand)




            flux_x = w_x_rho  * rho_i  + &
!                     rhow_x_rhou * rhou_i + &
!                     rhow_x_rhov * rhov_i + &
!                     rhow_x_rhow * rhow_i + &
                     w_x_p * p_i
            flux_y = w_y_rho  * rho_i  + &
!                     rhow_y_rhou * rhou_i + &
!                     rhow_y_rhov * rhov_i + &
!                     rhow_y_rhow * rhow_i + &
                     w_y_p * p_i
            flux_z = w_z_rho  * rho_i  + &
!                     rhow_z_rhou * rhou_i + &
!                     rhow_z_rhov * rhov_i + &
!                     rhow_z_rhow * rhow_i + &
                     w_z_p * p_i

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,iw_i,integrand)





            !=================================================
            ! Energy flux
            !=================================================

            flux_x = p_x_rho  * rho_r  + &
!                     rhoE_x_rhou * rhou_r + &
!                     rhoE_x_rhov * rhov_r + &
!                     rhoE_x_rhow * rhow_r + &
                     p_x_p * p_r
            flux_y = p_y_rho  * rho_r  + &
!                     rhoE_y_rhou * rhou_r + &
!                     rhoE_y_rhov * rhov_r + &
!                     rhoE_y_rhow * rhow_r + &
                     p_y_p * p_r
            flux_z = p_z_rho  * rho_r  + &
!                     rhoE_z_rhou * rhou_r + &
!                     rhoE_z_rhov * rhov_r + &
!                     rhoE_z_rhow * rhow_r + &
                     p_z_p * p_r

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,ip_r,integrand)


            flux_x = p_x_rho  * rho_i  + &
!                     rhoE_x_rhou * rhou_i + &
!                     rhoE_x_rhov * rhov_i + &
!                     rhoE_x_rhow * rhow_i + &
                     p_x_p * p_i
            flux_y = p_y_rho  * rho_i  + &
!                     rhoE_y_rhou * rhou_i + &
!                     rhoE_y_rhov * rhov_i + &
!                     rhoE_y_rhow * rhow_i + &
                     p_y_p * p_i
            flux_z = p_z_rho  * rho_i  + &
!                     rhoE_z_rhou * rhou_i + &
!                     rhoE_z_rhov * rhov_i + &
!                     rhoE_z_rhow * rhow_i + &
                     p_z_p * p_i

            integrand = flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3)

            call integrate_boundary_scalar_flux(mesh,sdata,face,fcn,ip_i,integrand)




            end associate

        end associate

    end subroutine compute
    !*********************************************************************************************************






end module bc_primlineuler_wall
