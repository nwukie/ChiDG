module PRIMLINEULER_boundary_average_advective_flux_imag
    use mod_kinds,              only: rk,ik
    use mod_constants,          only:  ZERO, ONE, TWO, HALF, ME, NEIGHBOR

    use type_boundary_flux,     only: boundary_flux_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D

    use PRIMLINEULER_properties,    only: PRIMLINEULER_properties_t
    use mod_primitive_linearized_euler
    implicit none
    private




    !> Implementation of the Euler boundary average flux
    !!
    !!  - At a boundary interface, the solution states Q- and Q+ exists on opposite 
    !!    sides of the boundary. The average flux is computed as Favg = 1/2(F(Q-) + F(Q+))
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/15/2016
    !!
    !------------------------------------------------------------------------------------------
    type, extends(boundary_flux_t), public :: PRIMLINEULER_boundary_average_advective_flux_imag_t

    contains

        procedure  :: compute

    end type PRIMLINEULER_boundary_average_advective_flux_imag_t
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
    subroutine compute(self,worker,prop)
        class(PRIMLINEULER_boundary_average_advective_flux_imag_t), intent(in)      :: self
        type(chidg_worker_t),                                       intent(inout)   :: worker
        class(properties_t),                                        intent(inout)   :: prop


        ! Equation indices
        integer(ik)     :: irho
        integer(ik)     :: iu
        integer(ik)     :: iv
        integer(ik)     :: iw
        integer(ik)     :: ip


        integer(ik)     :: ifcn, idonor, iblk, igq


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            rho_m,   rho_p,                         &
            u_m,     u_p,                           &
            v_m,     v_p,                           &
            w_m,     w_p,                           &
            p_m,     p_p,                           &
            flux_x_m,   flux_y_m,   flux_z_m,       &
            flux_x_p,   flux_y_p,   flux_z_p,       &
            flux_x,     flux_y,     flux_z,         &
            integrand

        real(rk),   allocatable, dimension(:)   ::  &
            normx, normy, normz

        irho = prop%get_eqn_index("rho_i")
        iu   = prop%get_eqn_index("u_i")
        iv   = prop%get_eqn_index("v_i")
        iw   = prop%get_eqn_index("w_i")
        ip   = prop%get_eqn_index("p_i")



        !
        ! Interpolate solution to quadrature nodes
        !
        rho_m = worker%interpolate(irho, 'value', ME)
        rho_p = worker%interpolate(irho, 'value', NEIGHBOR)

        u_m   = worker%interpolate(iu, 'value', ME)
        u_p   = worker%interpolate(iu, 'value', NEIGHBOR)

        v_m   = worker%interpolate(iv, 'value', ME)
        v_p   = worker%interpolate(iv, 'value', NEIGHBOR)

        w_m   = worker%interpolate(iw, 'value', ME)
        w_p   = worker%interpolate(iw, 'value', NEIGHBOR)

        p_m   = worker%interpolate(ip, 'value', ME)
        p_p   = worker%interpolate(ip, 'value', NEIGHBOR)




        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)




        !================================
        !       MASS FLUX
        !================================
        flux_x_m = rho_x_rho  * rho_m  + &
                   rho_x_u * u_m + &
                   rho_x_v * v_m + &
                   rho_x_w * w_m + &
                   rho_x_p * p_m
        flux_y_m = rho_y_rho  * rho_m  + &
                   rho_y_u * u_m + &
                   rho_y_v * v_m + &
                   rho_y_w * w_m + &
                   rho_y_p * p_m 
        flux_z_m = rho_z_rho  * rho_m  + &
                   rho_z_u * u_m + &
                   rho_z_v * v_m + &
                   rho_z_w * w_m + &
                   rho_z_p * p_m 


        flux_x_p = rho_x_rho  * rho_p  + &
                   rho_x_u * u_p + &
                   rho_x_v * v_p + &
                   rho_x_w * w_p + &
                   rho_x_p * p_p
        flux_y_p = rho_y_rho  * rho_p  + &
                   rho_y_u * u_p + &
                   rho_y_v * v_p + &
                   rho_y_w * w_p + &
                   rho_y_p * p_p 
        flux_z_p = rho_z_rho  * rho_p  + &
                   rho_z_u * u_p + &
                   rho_z_v * v_p + &
                   rho_z_w * w_p + &
                   rho_z_p * p_p 



        flux_x = (flux_x_m + flux_x_p)
        flux_y = (flux_y_m + flux_y_p)
        flux_z = (flux_z_m + flux_z_p)


        ! dot with normal vector
        integrand = HALF*(flux_x*normx + flux_y*normy + flux_z*normz)

        call worker%integrate_boundary(irho,integrand)


        !================================
        !       X-MOMENTUM FLUX
        !================================
        flux_x_m = u_x_rho  * rho_m  + &
                   u_x_u * u_m + &
                   u_x_v * v_m + &
                   u_x_w * w_m + &
                   u_x_p * p_m
        flux_y_m = u_y_rho  * rho_m  + &
                   u_y_u * u_m + &
                   u_y_v * v_m + &
                   u_y_w * w_m + &
                   u_y_p * p_m 
        flux_z_m = u_z_rho  * rho_m  + &
                   u_z_u * u_m + &
                   u_z_v * v_m + &
                   u_z_w * w_m + &
                   u_z_p * p_m 

        flux_x_p = u_x_rho  * rho_p  + &
                   u_x_u * u_p + &
                   u_x_v * v_p + &
                   u_x_w * w_p + &
                   u_x_p * p_p
        flux_y_p = u_y_rho  * rho_p  + &
                   u_y_u * u_p + &
                   u_y_v * v_p + &
                   u_y_w * w_p + &
                   u_y_p * p_p 
        flux_z_p = u_z_rho  * rho_p  + &
                   u_z_u * u_p + &
                   u_z_v * v_p + &
                   u_z_w * w_p + &
                   u_z_p * p_p 



        flux_x = (flux_x_m + flux_x_p)
        flux_y = (flux_y_m + flux_y_p)
        flux_z = (flux_z_m + flux_z_p)


        ! dot with normal vector
        integrand = HALF*(flux_x*normx + flux_y*normy + flux_z*normz)

        call worker%integrate_boundary(iu,integrand)


        !================================
        !       Y-MOMENTUM FLUX
        !================================
        flux_x_m = v_x_rho  * rho_m  + &
                   v_x_u * u_m + &
                   v_x_v * v_m + &
                   v_x_w * w_m + &
                   v_x_p * p_m
        flux_y_m = v_y_rho  * rho_m  + &
                   v_y_u * u_m + &
                   v_y_v * v_m + &
                   v_y_w * w_m + &
                   v_y_p * p_m 
        flux_z_m = v_z_rho  * rho_m  + &
                   v_z_u * u_m + &
                   v_z_v * v_m + &
                   v_z_w * w_m + &
                   v_z_p * p_m 


        flux_x_p = v_x_rho  * rho_p  + &
                   v_x_u * u_p + &
                   v_x_v * v_p + &
                   v_x_w * w_p + &
                   v_x_p * p_p
        flux_y_p = v_y_rho  * rho_p  + &
                   v_y_u * u_p + &
                   v_y_v * v_p + &
                   v_y_w * w_p + &
                   v_y_p * p_p 
        flux_z_p = v_z_rho  * rho_p  + &
                   v_z_u * u_p + &
                   v_z_v * v_p + &
                   v_z_w * w_p + &
                   v_z_p * p_p 



        flux_x = (flux_x_m + flux_x_p)
        flux_y = (flux_y_m + flux_y_p)
        flux_z = (flux_z_m + flux_z_p)


        ! dot with normal vector
        integrand = HALF*(flux_x*normx + flux_y*normy + flux_z*normz)

        call worker%integrate_boundary(iv,integrand)

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
        flux_z_m = w_z_rho  * rho_m  + &
                   w_z_u    * u_m + &
                   w_z_v    * v_m + &
                   w_z_w    * w_m + &
                   w_z_p    * p_m 


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
        flux_z_p = w_z_rho  * rho_p  + &
                   w_z_u    * u_p + &
                   w_z_v    * v_p + &
                   w_z_w    * w_p + &
                   w_z_p    * p_p 


        flux_x = (flux_x_m + flux_x_p)
        flux_y = (flux_y_m + flux_y_p)
        flux_z = (flux_z_m + flux_z_p)


        ! dot with normal vector
        integrand = HALF*(flux_x*normx + flux_y*normy + flux_z*normz)

        call worker%integrate_boundary(iw,integrand)




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
        flux_z_m = p_z_rho  * rho_m  + &
                   p_z_u    * u_m + &
                   p_z_v    * v_m + &
                   p_z_w    * w_m + &
                   p_z_p    * p_m 

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
        flux_z_p = p_z_rho  * rho_p  + &
                   p_z_u    * u_p + &
                   p_z_v    * v_p + &
                   p_z_w    * w_p + &
                   p_z_p    * p_p 


        flux_x = (flux_x_m + flux_x_p)
        flux_y = (flux_y_m + flux_y_p)
        flux_z = (flux_z_m + flux_z_p)


        ! dot with normal vector
        integrand = HALF*(flux_x*normx + flux_y*normy + flux_z*normz)

        call worker%integrate_boundary(ip,integrand)


    end subroutine compute
    !************************************************************************************************************












end module PRIMLINEULER_boundary_average_advective_flux_imag
