module PRIMLINEULER_volume_advective_flux_imag
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF,ZERO, &
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


    !>  Volume advective flux for Linearized Euler equations - imaginary.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!
    !----------------------------------------------------------------------------------
    type, extends(volume_flux_t), public :: PRIMLINEULER_volume_advective_flux_imag_t


    contains

        procedure  :: compute

    end type PRIMLINEULER_volume_advective_flux_imag_t
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
        class(PRIMLINEULER_volume_advective_flux_imag_t),   intent(in)      :: self
        type(mesh_t),                                   intent(in)      :: mesh(:)
        type(solverdata_t),                             intent(inout)   :: sdata
        class(properties_t),                            intent(inout)   :: prop
        integer(ik),                                    intent(in)      :: idom, ielem, iblk


        ! Equation indices
        integer(ik)    :: irho
        integer(ik)    :: iu
        integer(ik)    :: iv
        integer(ik)    :: iw
        integer(ik)    :: ip


        integer(ik)    :: iseed, i, idonor
        type(seed_t)   :: seed




        type(AD_D), dimension(mesh(idom)%elems(ielem)%gq%vol%nnodes)      ::  &
                    rho, u, v, w, p,                        &
                    flux_x, flux_y, flux_z


        idonor = 0


        !
        ! Get equation indices
        !
        irho  = prop%get_eqn_index("rho_i")
        iu = prop%get_eqn_index("u_i")
        iv = prop%get_eqn_index("v_i")
        iw = prop%get_eqn_index("w_i")
        ip = prop%get_eqn_index("p_i")





        !
        ! Get neighbor face and seed element for derivatives
        !
        seed = compute_seed(mesh,idom,ielem,iblk,idonor,iblk)




        !
        ! Interpolate solution to quadrature nodes
        !
        call interpolate_element(mesh,sdata%q,idom,ielem,irho, rho, seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,iu,u,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,iv,v,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,iw,w,seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,ip,p,seed)




        !===========================
        !        MASS FLUX
        !===========================
        flux_x = rho_x_rho  * rho  + &
                 rho_x_u    * u    + &
                 rho_x_v    * v    + &
                 rho_x_w    * w    + &
                 rho_x_p    * p
        flux_y = rho_y_rho  * rho  + &
                 rho_y_u    * u    + &
                 rho_y_v    * v    + &
                 rho_y_w    * w    + &
                 rho_y_p    * p
        flux_z = rho_z_rho  * rho  + &
                 rho_z_u    * u    + &
                 rho_z_v    * v    + &
                 rho_z_w    * w    + &
                 rho_z_p    * p

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irho,iblk,flux_x,flux_y,flux_z)


        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux_x = u_x_rho  * rho  + &
                 u_x_u    * u    + &
                 u_x_v    * v    + &
                 u_x_w    * w    + &
                 u_x_p    * p
        flux_y = u_y_rho  * rho  + &
                 u_y_u    * u    + &
                 u_y_v    * v    + &
                 u_y_w    * w    + &
                 u_y_p    * p
        flux_z = u_z_rho  * rho  + &
                 u_z_u    * u    + &
                 u_z_v    * v    + &
                 u_z_w    * w    + &
                 u_z_p    * p

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,iu,iblk,flux_x,flux_y,flux_z)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux_x = v_x_rho  * rho  + &
                 v_x_u    * u    + &
                 v_x_v    * v    + &
                 v_x_w    * w    + &
                 v_x_p    * p
        flux_y = v_y_rho  * rho  + &
                 v_y_u    * u    + &
                 v_y_v    * v    + &
                 v_y_w    * w    + &
                 v_y_p    * p
        flux_z = v_z_rho  * rho  + &
                 v_z_u    * u    + &
                 v_z_v    * v    + &
                 v_z_w    * w    + &
                 v_z_p    * p

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,iv,iblk,flux_x,flux_y,flux_z)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        flux_x = w_x_rho  * rho  + &
                 w_x_u    * u    + &
                 w_x_v    * v    + &
                 w_x_w    * w    + &
                 w_x_p    * p
        flux_y = w_y_rho  * rho  + &
                 w_y_u    * u    + &
                 w_y_v    * v    + &
                 w_y_w    * w    + &
                 w_y_p    * p
        flux_z = w_z_rho  * rho  + &
                 w_z_u    * u    + &
                 w_z_v    * v    + &
                 w_z_w    * w    + &
                 w_z_p    * p


        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,iw,iblk,flux_x,flux_y,flux_z)

        !============================
        !       ENERGY FLUX
        !============================
        flux_x = p_x_rho  * rho  + &
                 p_x_u    * u    + &
                 p_x_v    * v    + &
                 p_x_w    * w    + &
                 p_x_p    * p
        flux_y = p_y_rho  * rho  + &
                 p_y_u    * u    + &
                 p_y_v    * v    + &
                 p_y_w    * w    + &
                 p_y_p    * p
        flux_z = p_z_rho  * rho  + &
                 p_z_u    * u    + &
                 p_z_v    * v    + &
                 p_z_w    * w    + &
                 p_z_p    * p

        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,ip,iblk,flux_x,flux_y,flux_z)

    end subroutine compute
    !******************************************************************************************************






end module PRIMLINEULER_volume_advective_flux_imag
