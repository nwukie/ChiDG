module PRIMLINEULER_volume_advective_flux_real
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,THREE,FOUR,FIVE,EIGHT,NINE,HALF,ZERO,PI

    use type_volume_flux,       only: volume_flux_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
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
    subroutine compute(self,worker,prop)
        class(PRIMLINEULER_volume_advective_flux_real_t),   intent(in)      :: self
        type(chidg_worker_t),                               intent(inout)   :: worker
        class(properties_t),                                intent(inout)   :: prop


        ! Equation indices
        integer(ik)    :: irho_r, iu_r, iv_r, iw_r, ip_r
        integer(ik)    :: irho_i, iu_i, iv_i, iw_i, ip_i

        integer(ik) :: igq

        type(AD_D), allocatable, dimension(:)   ::  &
            rho_r, u_r, v_r, w_r, p_r,              &
            rho_i, u_i, v_i, w_i, p_i,              &
            flux_x, flux_y, flux_z

        real(rk),   allocatable, dimension(:)   ::  &
            x, y, z, r, sigma_x, sigma_y, sigma_z, fcn

        logical :: inA = .false.
        logical :: inB = .false.
        logical :: inC = .false.
        logical :: inD = .false.
        logical :: inE = .false.
        logical :: inF = .false.



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
        ! Interpolate solution to quadrature nodes
        !
        rho_r = worker%interpolate(irho_r, 'value')
        u_r   = worker%interpolate(iu_r,   'value')
        v_r   = worker%interpolate(iv_r,   'value')
        w_r   = worker%interpolate(iw_r,   'value')
        p_r   = worker%interpolate(ip_r,   'value')

        rho_i = worker%interpolate(irho_i, 'value')
        u_i   = worker%interpolate(iu_i,   'value')
        v_i   = worker%interpolate(iv_i,   'value')
        w_i   = worker%interpolate(iw_i,   'value')
        p_i   = worker%interpolate(ip_i,   'value')




        !
        ! Get coordinates
        !
        x = worker%x('boundary')
        y = worker%y('boundary')
        z = worker%z('boundary')
        r = sqrt(y**TWO + z**TWO)

        !
        ! Compute PML Layers
        !
        do igq = 1,size(x)


            ! Munt duct
            inA = ( x(igq) < -THREE + thickness ) .and. ( r(igq) > 1.212_rk )
            inB = ( x(igq) >  THREE - thickness )
            inC = ( y(igq) < -2.121_rk + thickness )
            inD = ( y(igq) >  2.121_rk - thickness )
            inE = ( z(igq) < -2.121_rk + thickness )
            inF = ( z(igq) >  2.121_rk - thickness )

!            ! Monopole
!            inA = ( x(igq) < -3._rk + thickness )
!            inB = ( x(igq) >  3._rk - thickness )
!            inC = ( y(igq) < -2.121_rk + thickness )
!            inD = ( y(igq) >  2.121_rk - thickness )
!            inE = ( z(igq) < -2.121_rk + thickness )
!            inF = ( z(igq) >  2.121_rk - thickness )


!            inA = .false.
!            inB = .false.
!            inC = .false.
!            inD = .false.
!            inE = .false.
!            inF = .false.


            ! X-PML
            if ( inA ) then
                fcn(igq)     =  abs( ( x(igq) - (-3._rk+thickness) ) / thickness )**TWO
                sigma_x(igq) = eps * fcn(igq)
            else if ( inB ) then
                fcn(igq)     =  abs( ( x(igq) - ( 3._rk-thickness) ) / thickness )**TWO
                sigma_x(igq) = eps * fcn(igq)
            else
                sigma_x(igq) = ZERO
            end if


            ! Y-PML
            if ( inC ) then
                fcn(igq)     =  abs( ( y(igq) - (-2.121_rk+thickness) ) / thickness )**TWO
                sigma_y(igq) = eps * fcn(igq)
            else if ( inD ) then
                fcn(igq)     =  abs( ( y(igq) - ( 2.121_rk-thickness) ) / thickness )**TWO
                sigma_y(igq) = eps * fcn(igq)
            else
                sigma_y(igq) = ZERO
            end if


            ! Z-PML
            if ( inE ) then
                fcn(igq)     =  abs( ( z(igq) - (-2.121_rk+thickness) ) / thickness )**TWO
                sigma_z(igq) = eps * fcn(igq)
            else if ( inF ) then
                fcn(igq)     =  abs( ( z(igq) - ( 2.121_rk-thickness) ) / thickness )**TWO
                sigma_z(igq) = eps * fcn(igq)
            else
                sigma_z(igq) = ZERO
            end if

        end do













        !===========================
        !        MASS FLUX
        !===========================
        flux_x = (rho_x_rho  * rho_r  + &
                  rho_x_u    * u_r    + &
                  rho_x_v    * v_r    + &
                  rho_x_w    * w_r    + &
                  rho_x_p    * p_r)*(ONE-sigma_y*sigma_z/(omega*omega))  - &
                 ! PML
                 (rho_x_rho * rho_i  + &
                  rho_x_u    * u_i    + &
                  rho_x_v    * v_i    + &
                  rho_x_w    * w_i    + &
                  rho_x_p    * p_i)*((sigma_y+sigma_z)/omega)

        flux_y = (rho_y_rho  * rho_r  + &
                  rho_y_u    * u_r    + &
                  rho_y_v    * v_r    + &
                  rho_y_w    * w_r    + &
                  rho_y_p    * p_r)*(ONE-sigma_x*sigma_z/(omega*omega))  - &
                 ! PML
                 (rho_y_rho * rho_i  + &
                  rho_y_u    * u_i    + &
                  rho_y_v    * v_i    + &
                  rho_y_w    * w_i    + &
                  rho_y_p    * p_i)*((sigma_x+sigma_z)/omega)

        flux_z = (rho_z_rho  * rho_r  + &
                  rho_z_u    * u_r    + &
                  rho_z_v    * v_r    + &
                  rho_z_w    * w_r    + &
                  rho_z_p    * p_r)*(ONE-sigma_x*sigma_y/(omega*omega))  - &
                 ! PML
                 (rho_z_rho * rho_i  + &
                  rho_z_u    * u_i    + &
                  rho_z_v    * v_i    + &
                  rho_z_w    * w_i    + &
                  rho_z_p    * p_i)*((sigma_x+sigma_y)/omega)

        call worker%integrate_volume(irho_r,flux_x,flux_y,flux_z)



        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux_x = (u_x_rho  * rho_r  + &
                 u_x_u    * u_r    + &
                 u_x_v    * v_r    + &
                 u_x_w    * w_r    + &
                 u_x_p    * p_r)*(ONE-sigma_y*sigma_z/(omega*omega))  - &
                 ! PML
                 (u_x_rho  * rho_i  + &
                  u_x_u    * u_i    + &
                  u_x_v    * v_i    + &
                  u_x_w    * w_i    + &
                  u_x_p    * p_i)*((sigma_y+sigma_z)/omega)

        flux_y = (u_y_rho  * rho_r  + &
                 u_y_u    * u_r    + &
                 u_y_v    * v_r    + &
                 u_y_w    * w_r    + &
                 u_y_p    * p_r)*(ONE-sigma_x*sigma_z/(omega*omega))  - &
                 ! PML
                 (u_y_rho  * rho_i  + &
                  u_y_u    * u_i    + &
                  u_y_v    * v_i    + &
                  u_y_w    * w_i    + &
                  u_y_p    * p_i)*((sigma_x+sigma_z)/omega)

        flux_z = (u_z_rho  * rho_r  + &
                 u_z_u    * u_r    + &
                 u_z_v    * v_r    + &
                 u_z_w    * w_r    + &
                 u_z_p    * p_r)*(ONE-sigma_x*sigma_y/(omega*omega))  - &
                 ! PML
                 (u_z_rho * rho_i  + &
                  u_z_u    * u_i    + &
                  u_z_v    * v_i    + &
                  u_z_w    * w_i    + &
                  u_z_p    * p_i)*((sigma_x+sigma_y)/omega)

        call worker%integrate_volume(iu_r,flux_x,flux_y,flux_z)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux_x = (v_x_rho  * rho_r  + &
                 v_x_u    * u_r    + &
                 v_x_v    * v_r    + &
                 v_x_w    * w_r    + &
                 v_x_p    * p_r)*(ONE-sigma_y*sigma_z/(omega*omega))  - &
                 ! PML
                 (v_x_rho  * rho_i  + &
                  v_x_u    * u_i    + &
                  v_x_v    * v_i    + &
                  v_x_w    * w_i    + &
                  v_x_p    * p_i)*((sigma_y+sigma_z)/omega)

        flux_y = (v_y_rho  * rho_r  + &
                 v_y_u    * u_r    + &
                 v_y_v    * v_r    + &
                 v_y_w    * w_r    + &
                 v_y_p    * p_r)*(ONE-sigma_x*sigma_z/(omega*omega))  - &
                 ! PML
                 (v_y_rho * rho_i  + &
                  v_y_u    * u_i    + &
                  v_y_v    * v_i    + &
                  v_y_w    * w_i    + &
                  v_y_p    * p_i)*((sigma_x+sigma_z)/omega)

        flux_z = (v_z_rho  * rho_r  + &
                 v_z_u    * u_r    + &
                 v_z_v    * v_r    + &
                 v_z_w    * w_r    + &
                 v_z_p    * p_r)*(ONE-sigma_x*sigma_y/(omega*omega))  - &
                 ! PML
                 (v_z_rho * rho_i  + &
                  v_z_u    * u_i    + &
                  v_z_v    * v_i    + &
                  v_z_w    * w_i    + &
                  v_z_p    * p_i)*((sigma_x+sigma_y)/omega)

        call worker%integrate_volume(iv_r,flux_x,flux_y,flux_z)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        flux_x = (w_x_rho  * rho_r  + &
                 w_x_u    * u_r    + &
                 w_x_v    * v_r    + &
                 w_x_w    * w_r    + &
                 w_x_p    * p_r)*(ONE-sigma_y*sigma_z/(omega*omega))  - &
                 ! PML
                 (w_x_rho * rho_i  + &
                  w_x_u   * u_i    + &
                  w_x_v   * v_i    + &
                  w_x_w   * w_i    + &
                  w_x_p   * p_i)*((sigma_y+sigma_z)/omega)

        flux_y = (w_y_rho  * rho_r  + &
                 w_y_u    * u_r    + &
                 w_y_v    * v_r    + &
                 w_y_w    * w_r    + &
                 w_y_p    * p_r)*(ONE-sigma_x*sigma_z/(omega*omega))  - &
                 ! PML
                 (w_y_rho * rho_i  + &
                  w_y_u    * u_i    + &
                  w_y_v    * v_i    + &
                  w_y_w    * w_i    + &
                  w_y_p    * p_i)*((sigma_x+sigma_z)/omega)

        flux_z = (w_z_rho  * rho_r  + &
                 w_z_u    * u_r    + &
                 w_z_v    * v_r    + &
                 w_z_w    * w_r    + &
                 w_z_p    * p_r)*(ONE-sigma_x*sigma_y/(omega*omega))  - &
                 ! PML
                 (w_z_rho * rho_i  + &
                  w_z_u    * u_i    + &
                  w_z_v    * v_i    + &
                  w_z_w    * w_i    + &
                  w_z_p    * p_i)*((sigma_x+sigma_y)/omega)

        call worker%integrate_volume(iw_r,flux_x,flux_y,flux_z)

        !============================
        !       ENERGY FLUX
        !============================
        flux_x = (p_x_rho  * rho_r  + &
                 p_x_u    * u_r    + &
                 p_x_v    * v_r    + &
                 p_x_w    * w_r    + &
                 p_x_p    * p_r)*(ONE-sigma_y*sigma_z/(omega*omega))  - &
                 ! PML
                 (p_x_rho  * rho_i  + &
                  p_x_u    * u_i    + &
                  p_x_v    * v_i    + &
                  p_x_w    * w_i    + &
                  p_x_p    * p_i)*((sigma_y+sigma_z)/omega)

        flux_y = (p_y_rho  * rho_r  + &
                 p_y_u    * u_r    + &
                 p_y_v    * v_r    + &
                 p_y_w    * w_r    + &
                 p_y_p    * p_r)*(ONE-sigma_x*sigma_z/(omega*omega))  - &
                 ! PML
                 (p_y_rho * rho_i  + &
                  p_y_u    * u_i    + &
                  p_y_v    * v_i    + &
                  p_y_w    * w_i    + &
                  p_y_p    * p_i)*((sigma_x+sigma_z)/omega)

        flux_z = (p_z_rho  * rho_r  + &
                 p_z_u    * u_r    + &
                 p_z_v    * v_r    + &
                 p_z_w    * w_r    + &
                 p_z_p    * p_r)*(ONE-sigma_x*sigma_y/(omega*omega))  - &
                 ! PML
                 (p_z_rho * rho_i  + &
                  p_z_u    * u_i    + &
                  p_z_v    * v_i    + &
                  p_z_w    * w_i    + &
                  p_z_p    * p_i)*((sigma_x+sigma_y)/omega)

        call worker%integrate_volume(ip_r,flux_x,flux_y,flux_z)

    end subroutine compute
    !******************************************************************************************************






end module PRIMLINEULER_volume_advective_flux_real
