module PRIMLINEULER_volume_advective_flux_real
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,THREE,HALF,ZERO,PI,&
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use type_mesh,              only: mesh_t
    use atype_volume_flux,      only: volume_flux_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    
    use mod_interpolate,        only: interpolate_element
    use mod_integrate,          only: integrate_volume_flux
    use mod_DNAD_tools
    use DNAD_D

    use PRIMLINEULER_properties,    only: PRIMLINEULER_properties_t
    use mod_primitive_linearized_euler
    implicit none

    private


    !>  Volume advective flux for Linearized Euler equations - real.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!
    !----------------------------------------------------------------------------------
    type, extends(volume_flux_t), public :: PRIMLINEULER_volume_advective_flux_real_t


    contains

        procedure  :: compute

    end type PRIMLINEULER_volume_advective_flux_real_t
    !***********************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,idom,ielem,iblk)
        class(PRIMLINEULER_volume_advective_flux_real_t),   intent(in)      :: self
        type(mesh_t),                                   intent(in)      :: mesh(:)
        type(solverdata_t),                             intent(inout)   :: sdata
        class(properties_t),                            intent(inout)   :: prop
        integer(ik),                                    intent(in)      :: idom, ielem, iblk

        ! Equation indices
        integer(ik)    :: irho_r, iu_r, iv_r, iw_r, ip_r
        integer(ik)    :: irho_i, iu_i, iv_i, iw_i, ip_i


        integer(ik)    :: iseed, i, idonor, igq
        type(seed_t)   :: seed

        real(rk)    :: gam, omega, thickness, eps, kappa

        type(AD_D), dimension(mesh(idom)%elems(ielem)%gq%vol%nnodes)      ::  &
                    rho_r, u_r, v_r, w_r, p_r,                        &
                    rho_i, u_i, v_i, w_i, p_i,                        &
                    flux_x, flux_y, flux_z

        real(rk), dimension(mesh(idom)%elems(ielem)%gq%vol%nnodes)      ::  &
                    x, y, sigma_x, sigma_y, sigma, fcn

        logical :: inA = .false.
        logical :: inB = .false.
        logical :: inC = .false.




        idonor = 0

        !
        ! Get equation indices
        !
        irho_i = prop%get_eqn_index("rho_i")
        iu_i   = prop%get_eqn_index("u_i")
        iv_i   = prop%get_eqn_index("v_i")
        iw_i   = prop%get_eqn_index("w_i")
        ip_i   = prop%get_eqn_index("p_i")

        irho_r = prop%get_eqn_index("rho_r")
        iu_r   = prop%get_eqn_index("u_r")
        iv_r   = prop%get_eqn_index("v_r")
        iw_r   = prop%get_eqn_index("w_r")
        ip_r   = prop%get_eqn_index("p_r")





        !
        ! Get neighbor face and seed element for derivatives
        !
        seed = compute_seed(mesh,idom,ielem,iblk,idonor,iblk)

        !
        ! Interpolate solution to quadrature nodes
        !
        call interpolate_element(mesh,sdata%q,idom,ielem,irho_r, rho_r, seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,iu_r,u_r,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,iv_r,v_r,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,iw_r,w_r,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,ip_r,p_r,seed)


        call interpolate_element(mesh,sdata%q,idom,ielem,irho_i, rho_i, seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,iu_i,u_i,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,iv_i,v_i,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,iw_i,w_i,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,ip_i,p_i,seed)



        !
        ! Gamma
        !
        gam = 1.4_rk
        omega = 956._rk * TWO * PI



        !
        ! Absorbing layer
        !
        thickness = HALF
        eps       = 100._rk
        kappa     = 1._rk

        ! Get coordinates
        x = mesh(idom)%elems(ielem)%quad_pts(:)%c1_
        y = mesh(idom)%elems(ielem)%quad_pts(:)%c2_

        do igq = 1,size(x)

            inA = ( x(igq) > -THREE ) .and. ( x(igq) < -THREE + thickness ) .and. ( y(igq) > 1.2_rk )
            inB = ( y(igq) > 4.6_rk - thickness )  .and.  ( y(igq) < 4.6_rk )
            inC = ( x(igq) > 6.2_rk - thickness )

            if ( inA ) then
                fcn     = -(ONE/thickness)*x  +  (ONE - THREE/thickness)
!                sigma_x = eps*(ONE-exp(kappa*fcn**TWO))/(ONE-exp(kappa))
                sigma_x = fcn
            
            else if ( inC ) then
                fcn     =  (ONE/thickness)*x  +  (ONE - 6.2_rk/thickness)
!                sigma_x = eps*(ONE-exp(kappa*fcn**TWO))/(ONE-exp(kappa))
                sigma_x = fcn

            else
                sigma_x = ZERO

            end if


            if ( inB ) then
                fcn     =  (ONE/thickness)*y  +  (ONE - 4.6_rk/thickness)
!                sigma_y = eps*(ONE-exp(kappa*fcn**TWO))/(ONE-exp(kappa))
                sigma_y = fcn

            else
                sigma_y = ZERO
            end if

            sigma = sigma_x * sigma_y

        end do
















        !===========================
        !        MASS FLUX
        !===========================
        flux_x = rho_x_rho  * rho_r  + &
                 rho_x_u    * u_r    + &
                 rho_x_v    * v_r    + &
                 rho_x_w    * w_r    + &
                 rho_x_p    * p_r    - &
                 ! PML
                 (rho_x_rho * rho_i  + &
                 rho_x_u    * u_i    + &
                 rho_x_v    * v_i    + &
                 rho_x_w    * w_i    + &
                 rho_x_p    * p_i)*sigma_y/omega

        flux_y = rho_y_rho  * rho_r  + &
                 rho_y_u    * u_r    + &
                 rho_y_v    * v_r    + &
                 rho_y_w    * w_r    + &
                 rho_y_p    * p_r    - &
                 ! PML
                 (rho_y_rho * rho_i  + &
                 rho_y_u    * u_i    + &
                 rho_y_v    * v_i    + &
                 rho_y_w    * w_i    + &
                 rho_y_p    * p_i)*sigma_x/omega

        flux_z = rho_z_rho  * rho_r  + &
                 rho_z_u    * u_r    + &
                 rho_z_v    * v_r    + &
                 rho_z_w    * w_r    + &
                 rho_z_p    * p_r

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irho_r,iblk,flux_x,flux_y,flux_z)


        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux_x = u_x_rho  * rho_r  + &
                 u_x_u    * u_r    + &
                 u_x_v    * v_r    + &
                 u_x_w    * w_r    + &
                 u_x_p    * p_r    - &
                 ! PML
                 (u_x_rho * rho_i  + &
                 u_x_u    * u_i    + &
                 u_x_v    * v_i    + &
                 u_x_w    * w_i    + &
                 u_x_p    * p_i)*sigma_y/omega

        flux_y = u_y_rho  * rho_r  + &
                 u_y_u    * u_r    + &
                 u_y_v    * v_r    + &
                 u_y_w    * w_r    + &
                 u_y_p    * p_r    - &
                 ! PML
                 (u_y_rho * rho_i  + &
                 u_y_u    * u_i    + &
                 u_y_v    * v_i    + &
                 u_y_w    * w_i    + &
                 u_y_p    * p_i)*sigma_x/omega

        flux_z = u_z_rho  * rho_r  + &
                 u_z_u    * u_r    + &
                 u_z_v    * v_r    + &
                 u_z_w    * w_r    + &
                 u_z_p    * p_r

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,iu_r,iblk,flux_x,flux_y,flux_z)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux_x = v_x_rho  * rho_r  + &
                 v_x_u    * u_r    + &
                 v_x_v    * v_r    + &
                 v_x_w    * w_r    + &
                 v_x_p    * p_r    - &
                 ! PML
                 (v_x_rho * rho_i  + &
                 v_x_u    * u_i    + &
                 v_x_v    * v_i    + &
                 v_x_w    * w_i    + &
                 v_x_p    * p_i)*sigma_y/omega

        flux_y = v_y_rho  * rho_r  + &
                 v_y_u    * u_r    + &
                 v_y_v    * v_r    + &
                 v_y_w    * w_r    + &
                 v_y_p    * p_r    - &
                 ! PML
                 (v_y_rho * rho_i  + &
                 v_y_u    * u_i    + &
                 v_y_v    * v_i    + &
                 v_y_w    * w_i    + &
                 v_y_p    * p_i)*sigma_x/omega

        flux_z = v_z_rho  * rho_r  + &
                 v_z_u    * u_r    + &
                 v_z_v    * v_r    + &
                 v_z_w    * w_r    + &
                 v_z_p    * p_r

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,iv_r,iblk,flux_x,flux_y,flux_z)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        flux_x = w_x_rho  * rho_r  + &
                 w_x_u    * u_r    + &
                 w_x_v    * v_r    + &
                 w_x_w    * w_r    + &
                 w_x_p    * p_r    - &
                 ! PML
                 (w_x_rho * rho_i  + &
                 w_x_u    * u_i    + &
                 w_x_v    * v_i    + &
                 w_x_w    * w_i    + &
                 w_x_p    * p_i)*sigma_y/omega

        flux_y = w_y_rho  * rho_r  + &
                 w_y_u    * u_r    + &
                 w_y_v    * v_r    + &
                 w_y_w    * w_r    + &
                 w_y_p    * p_r    - &
                 ! PML
                 (w_y_rho * rho_i  + &
                 w_y_u    * u_i    + &
                 w_y_v    * v_i    + &
                 w_y_w    * w_i    + &
                 w_y_p    * p_i)*sigma_x/omega

        flux_z = w_z_rho  * rho_r  + &
                 w_z_u    * u_r    + &
                 w_z_v    * v_r    + &
                 w_z_w    * w_r    + &
                 w_z_p    * p_r

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,iw_r,iblk,flux_x,flux_y,flux_z)

        !============================
        !       ENERGY FLUX
        !============================
        flux_x = p_x_rho  * rho_r  + &
                 p_x_u    * u_r    + &
                 p_x_v    * v_r    + &
                 p_x_w    * w_r    + &
                 p_x_p    * p_r    - &
                 ! PML
                 (p_x_rho * rho_i  + &
                 p_x_u    * u_i    + &
                 p_x_v    * v_i    + &
                 p_x_w    * w_i    + &
                 p_x_p    * p_i)*sigma_y/omega

        flux_y = p_y_rho  * rho_r  + &
                 p_y_u    * u_r    + &
                 p_y_v    * v_r    + &
                 p_y_w    * w_r    + &
                 p_y_p    * p_r    - &
                 ! PML
                 (p_y_rho * rho_i  + &
                 p_y_u    * u_i    + &
                 p_y_v    * v_i    + &
                 p_y_w    * w_i    + &
                 p_y_p    * p_i)*sigma_x/omega

        flux_z = p_z_rho  * rho_r  + &
                 p_z_u    * u_r    + &
                 p_z_v    * v_r    + &
                 p_z_w    * w_r    + &
                 p_z_p    * p_r

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,ip_r,iblk,flux_x,flux_y,flux_z)

    end subroutine compute
    !******************************************************************************************************






end module PRIMLINEULER_volume_advective_flux_real
