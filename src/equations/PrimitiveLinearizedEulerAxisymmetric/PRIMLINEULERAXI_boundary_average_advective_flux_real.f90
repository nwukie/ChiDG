module PRIMLINEULERAXI_boundary_average_advective_flux_real
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES, ZERO, ONE, TWO, HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX, &
                                      LOCAL, NEIGHBOR

    use atype_boundary_flux,    only: boundary_flux_t
    use type_mesh,              only: mesh_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    use type_seed,              only: seed_t
    use type_face_info,         only: face_info_t
    use type_function_info,     only: function_info_t

    use mod_interpolate,        only: interpolate_face
    use mod_integrate,          only: integrate_boundary_scalar_flux
    use mod_DNAD_tools
    use DNAD_D

    use PRIMLINEULERAXI_properties,    only: PRIMLINEULERAXI_properties_t
    use mod_primitive_linearized_euler_axisymmetric
    implicit none



    private



    !> Implementation of the Euler boundary average flux
    !!
    !!  - At a boundary interface, the solution states Q- and Q+ exists on opposite 
    !!    sides of the boundary. The average flux is computed as Favg = 1/2(F(Q-) + F(Q+))
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !------------------------------------------------------------------------------------------
    type, extends(boundary_flux_t), public :: PRIMLINEULERAXI_boundary_average_advective_flux_real_t

    contains
        procedure  :: compute

    end type PRIMLINEULERAXI_boundary_average_advective_flux_real_t
    !*******************************************************************************************










contains



    !>   Boundary Flux routine for Euler
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/16/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,face_info,function_info)
        class(PRIMLINEULERAXI_boundary_average_advective_flux_real_t), intent(in)      :: self
        type(mesh_t),                                           intent(in)      :: mesh(:)
        type(solverdata_t),                                     intent(inout)   :: sdata
        class(properties_t),                                    intent(inout)   :: prop
        type(face_info_t),                                      intent(in)      :: face_info
        type(function_info_t),                                  intent(in)      :: function_info

        ! Equation indices
        integer(ik)     :: irho
        integer(ik)     :: iu
        integer(ik)     :: iv
        integer(ik)     :: iw
        integer(ik)     :: ip


        integer(ik) :: idom, ielem, iface
        integer(ik) :: ifcn, idonor, iblk, igq

        ! Storage at quadrature nodes
        type(AD_D), dimension(mesh(face_info%idomain)%faces(face_info%ielement,face_info%iface)%gq%face%nnodes)    :: &
                        rho_m,   rho_p,                                  &
                        u_m,     u_p,                                 &
                        v_m,     v_p,                                 &
                        w_m,     w_p,                                 &
                        p_m,     p_p,                                 &
                        flux_x_m,   flux_y_m,   flux_z_m,                   &
                        flux_x_p,   flux_y_p,   flux_z_p,                   &
                        flux_x,     flux_y,     flux_z,                     &
                        integrand

!        real(rk), dimension(mesh(face_info%idomain)%faces(face_info%ielement,face_info%iface)%gq%face%nnodes)    :: &
!                        y


        !===========================================================================
        ! NOTE: var_m signifies "minus" and would indicate a local element variable
        !       var_p signifies "plus"  and would indicate a neighbor element variable
        !===========================================================================
        irho  = prop%get_eqn_index("rho_r")
        iu = prop%get_eqn_index("u_r")
        iv = prop%get_eqn_index("v_r")
        iw = prop%get_eqn_index("w_r")
        ip = prop%get_eqn_index("p_r")

        idom  = face_info%idomain
        ielem = face_info%ielement
        iface = face_info%iface


        associate (norms => mesh(idom)%faces(ielem,iface)%norm, faces => mesh(idom)%faces, q => sdata%q)


            !
            ! Interpolate solution to quadrature nodes
            !
            call interpolate_face(mesh,face_info,q, irho,  rho_m,  LOCAL)
            call interpolate_face(mesh,face_info,q, irho,  rho_p,  NEIGHBOR)

            call interpolate_face(mesh,face_info,q, iu, u_m, LOCAL)
            call interpolate_face(mesh,face_info,q, iu, u_p, NEIGHBOR)

            call interpolate_face(mesh,face_info,q, iv, v_m, LOCAL)
            call interpolate_face(mesh,face_info,q, iv, v_p, NEIGHBOR)

            call interpolate_face(mesh,face_info,q, iw, w_m, LOCAL)
            call interpolate_face(mesh,face_info,q, iw, w_p, NEIGHBOR)

            call interpolate_face(mesh,face_info,q, ip, p_m, LOCAL)
            call interpolate_face(mesh,face_info,q, ip, p_p, NEIGHBOR)


!            !
!            ! Get y-coordinate
!            !
!            y = mesh(idom)%faces(ielem,iface)%quad_pts(:)%c2_



            !================================
            !       MASS FLUX
            !================================
            flux_x_m = rho_x_rho  * rho_m  + &
                       rho_x_u    * u_m + &
                       rho_x_v    * v_m + &
                       rho_x_w    * w_m + &
                       rho_x_p    * p_m
            flux_y_m = rho_y_rho  * rho_m  + &
                       rho_y_u    * u_m + &
                       rho_y_v    * v_m + &
                       rho_y_w    * w_m + &
                       rho_y_p    * p_m 
!            flux_z_m = rho_z_rho  * rho_m  + &
!                       rho_z_u    * u_m + &
!                       rho_z_v    * v_m + &
!                       rho_z_w    * w_m + &
!                       rho_z_p    * p_m 

            flux_x_p = rho_x_rho  * rho_p  + &
                       rho_x_u    * u_p + &
                       rho_x_v    * v_p + &
                       rho_x_w    * w_p + &
                       rho_x_p    * p_p
            flux_y_p = rho_y_rho  * rho_p  + &
                       rho_y_u    * u_p + &
                       rho_y_v    * v_p + &
                       rho_y_w    * w_p + &
                       rho_y_p    * p_p 
!            flux_z_p = rho_z_rho  * rho_p  + &
!                       rho_z_u    * u_p + &
!                       rho_z_v    * v_p + &
!                       rho_z_w    * w_p + &
!                       rho_z_p    * p_p 


            flux_x = (flux_x_m + flux_x_p)
            flux_y = (flux_y_m + flux_y_p)
            flux_z = flux_y
            flux_z = ZERO
!            flux_z = (flux_z_m + flux_z_p)



            ! dot with normal vector
            integrand = HALF*(flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,irho,integrand)




            !================================
            !       X-MOMENTUM FLUX
            !================================
            flux_x_m = u_x_rho  * rho_m  + &
                       u_x_u    * u_m + &
                       u_x_v    * v_m + &
                       u_x_w    * w_m + &
                       u_x_p    * p_m
            flux_y_m = u_y_rho  * rho_m  + &
                       u_y_u    * u_m + &
                       u_y_v    * v_m + &
                       u_y_w    * w_m + &
                       u_y_p    * p_m 
!            flux_z_m = u_z_rho  * rho_m  + &
!                       u_z_u    * u_m + &
!                       u_z_v    * v_m + &
!                       u_z_w    * w_m + &
!                       u_z_p    * p_m 

            flux_x_p = u_x_rho  * rho_p  + &
                       u_x_u    * u_p + &
                       u_x_v    * v_p + &
                       u_x_w    * w_p + &
                       u_x_p    * p_p
            flux_y_p = u_y_rho  * rho_p  + &
                       u_y_u    * u_p + &
                       u_y_v    * v_p + &
                       u_y_w    * w_p + &
                       u_y_p    * p_p 
!            flux_z_p = u_z_rho  * rho_p  + &
!                       u_z_u    * u_p + &
!                       u_z_v    * v_p + &
!                       u_z_w    * w_p + &
!                       u_z_p    * p_p 


            flux_x = (flux_x_m + flux_x_p)
            flux_y = (flux_y_m + flux_y_p)
            flux_z = flux_y
            flux_z = ZERO
!            flux_z = (flux_z_m + flux_z_p)


            ! dot with normal vector
            integrand = HALF*(flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,iu,integrand)



            !================================
            !       Y-MOMENTUM FLUX
            !================================
            flux_x_m = v_x_rho  * rho_m  + &
                       v_x_u    * u_m + &
                       v_x_v    * v_m + &
                       v_x_w    * w_m + &
                       v_x_p    * p_m
            flux_y_m = v_y_rho  * rho_m  + &
                       v_y_u    * u_m + &
                       v_y_v    * v_m + &
                       v_y_w    * w_m + &
                       v_y_p    * p_m 
!            flux_z_m = v_z_rho  * rho_m  + &
!                       v_z_u    * u_m + &
!                       v_z_v    * v_m + &
!                       v_z_w    * w_m + &
!                       v_z_p    * p_m 


            flux_x_p = v_x_rho  * rho_p  + &
                       v_x_u    * u_p + &
                       v_x_v    * v_p + &
                       v_x_w    * w_p + &
                       v_x_p    * p_p
            flux_y_p = v_y_rho  * rho_p  + &
                       v_y_u    * u_p + &
                       v_y_v    * v_p + &
                       v_y_w    * w_p + &
                       v_y_p    * p_p 
!            flux_z_p = v_z_rho  * rho_p  + &
!                       v_z_u    * u_p + &
!                       v_z_v    * v_p + &
!                       v_z_w    * w_p + &
!                       v_z_p    * p_p 



            flux_x = (flux_x_m + flux_x_p)
            flux_y = (flux_y_m + flux_y_p)
            flux_z = flux_y
            flux_z = ZERO
!            flux_z = (flux_z_m + flux_z_p)



            ! dot with normal vector
            integrand = HALF*(flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,iv,integrand)




            !================================
            !       Z-MOMENTUM FLUX
            !================================
            flux_x_m = w_x_rho  * rho_m  + &
                       w_x_u    * u_m + &
                       w_x_v    * v_m + &
                       w_x_w    * w_m + &
                       w_x_p    * p_m
            flux_y_m = w_y_rho  * rho_m  + &
                       w_y_u    * u_m + &
                       w_y_v    * v_m + &
                       w_y_w    * w_m + &
                       w_y_p    * p_m 
!            flux_z_m = w_z_rho  * rho_m  + &
!                       w_z_u    * u_m + &
!                       w_z_v    * v_m + &
!                       w_z_w    * w_m + &
!                       w_z_p    * p_m 


            flux_x_p = w_x_rho  * rho_p  + &
                       w_x_u    * u_p + &
                       w_x_v    * v_p + &
                       w_x_w    * w_p + &
                       w_x_p    * p_p
            flux_y_p = w_y_rho  * rho_p  + &
                       w_y_u    * u_p + &
                       w_y_v    * v_p + &
                       w_y_w    * w_p + &
                       w_y_p    * p_p 
!            flux_z_p = w_z_rho  * rho_p  + &
!                       w_z_u    * u_p + &
!                       w_z_v    * v_p + &
!                       w_z_w    * w_p + &
!                       w_z_p    * p_p 



            flux_x = (flux_x_m + flux_x_p)
            flux_y = (flux_y_m + flux_y_p)
            flux_z = flux_y
            flux_z = ZERO
!            flux_z = (flux_z_m + flux_z_p)



            ! dot with normal vector
            integrand = HALF*(flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,iw,integrand)


            !================================
            !          ENERGY FLUX
            !================================
            flux_x_m = p_x_rho  * rho_m  + &
                       p_x_u    * u_m + &
                       p_x_v    * v_m + &
                       p_x_w    * w_m + &
                       p_x_p    * p_m
            flux_y_m = p_y_rho  * rho_m  + &
                       p_y_u    * u_m + &
                       p_y_v    * v_m + &
                       p_y_w    * w_m + &
                       p_y_p    * p_m 
!            flux_z_m = p_z_rho  * rho_m  + &
!                       p_z_u    * u_m + &
!                       p_z_v    * v_m + &
!                       p_z_w    * w_m + &
!                       p_z_p    * p_m 

            flux_x_p = p_x_rho  * rho_p  + &
                       p_x_u    * u_p + &
                       p_x_v    * v_p + &
                       p_x_w    * w_p + &
                       p_x_p    * p_p
            flux_y_p = p_y_rho  * rho_p  + &
                       p_y_u    * u_p + &
                       p_y_v    * v_p + &
                       p_y_w    * w_p + &
                       p_y_p    * p_p 
!            flux_z_p = p_z_rho  * rho_p  + &
!                       p_z_u    * u_p + &
!                       p_z_v    * v_p + &
!                       p_z_w    * w_p + &
!                       p_z_p    * p_p 



            flux_x = (flux_x_m + flux_x_p)
            flux_y = (flux_y_m + flux_y_p)
            flux_z = flux_y
            flux_z = ZERO
!            flux_z = (flux_z_m + flux_z_p)


            ! dot with normal vector
            integrand = HALF*(flux_x*norms(:,1) + flux_y*norms(:,2) + flux_z*norms(:,3))

            call integrate_boundary_scalar_flux(mesh,sdata,face_info,function_info,ip,integrand)




        end associate

    end subroutine compute
    !********************************************************************************************************












end module PRIMLINEULERAXI_boundary_average_advective_flux_real
