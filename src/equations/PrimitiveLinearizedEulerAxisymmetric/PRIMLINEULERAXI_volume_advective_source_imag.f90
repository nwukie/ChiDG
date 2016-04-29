module PRIMLINEULERAXI_volume_advective_source_imag
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,THREE,FOUR,FIVE,EIGHT,NINE,HALF,ZERO, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG,PI

    use type_mesh,              only: mesh_t
    use atype_volume_flux,      only: volume_flux_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    
    use mod_interpolate,        only: interpolate_element
    use mod_integrate,          only: integrate_volume_source
    use mod_DNAD_tools
    use DNAD_D

    use PRIMLINEULERAXI_properties,                 only: PRIMLINEULERAXI_properties_t
    !use mod_primitive_linearized_euler_axisymmetric, only: omega, gam, thickness, eps, rhobar, pbar, mod_m
    use mod_primitive_linearized_euler_axisymmetric
    implicit none

    private



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !-----------------------------------------------------------------------------------
    type, extends(volume_flux_t), public :: PRIMLINEULERAXI_volume_advective_source_imag_t


    contains

        procedure  :: compute

    end type PRIMLINEULERAXI_volume_advective_source_imag_t
    !************************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,idom,ielem,iblk)
        class(PRIMLINEULERAXI_volume_advective_source_imag_t), intent(in)      :: self
        type(mesh_t),                                   intent(in)      :: mesh(:)
        type(solverdata_t),                             intent(inout)   :: sdata
        class(properties_t),                            intent(inout)   :: prop
        integer(ik),                                    intent(in)      :: idom, ielem, iblk

        ! Equation indices
        integer(ik)    :: irho_r,  irho_i
        integer(ik)    :: iu_r, iu_i
        integer(ik)    :: iv_r, iv_i
        integer(ik)    :: iw_r, iw_i
        integer(ik)    :: ip_r, ip_i


        integer(ik)    :: iseed, i, idonor, igq
        type(seed_t)   :: seed

!        real(rk)    :: gam, thickness, eps, kappa



        type(AD_D), dimension(mesh(idom)%elems(ielem)%gq%vol%nnodes)      ::    &
                    rho_r, u_r, v_r, w_r, p_r,                      & 
                    rho_i, u_i, v_i, w_i, p_i,                      &
                    p, H,                                                       &
                    flux

        real(rk), dimension(mesh(idom)%elems(ielem)%gq%vol%nnodes)      ::  &
                    x, y, sigma_x, sigma_y, sigma, fcn

        logical :: inA = .false.
        logical :: inB = .false.
        logical :: inC = .false.
        logical :: inD = .false.


        idonor = 0


        !-------------------------------------------------------------
        irho_r  = prop%get_eqn_index("rho_r")
        iu_r = prop%get_eqn_index("u_r")
        iv_r = prop%get_eqn_index("v_r")
        iw_r = prop%get_eqn_index("w_r")
        ip_r = prop%get_eqn_index("p_r")

        irho_i  = prop%get_eqn_index("rho_i")
        iu_i = prop%get_eqn_index("u_i")
        iv_i = prop%get_eqn_index("v_i")
        iw_i = prop%get_eqn_index("w_i")
        ip_i = prop%get_eqn_index("p_i")



        !
        ! Gamma
        !
        !gam = 1.4_rk
        !omega = 956._rk * TWO * PI
        !omega = 1200._rk * TWO * PI





        !
        ! Absorbing layer
        !
        !thickness = 0.7_rk
        !eps       = 1700._rk
        !kappa     = 1._rk

        ! Get coordinates
        x = mesh(idom)%elems(ielem)%quad_pts(:)%c1_
        y = mesh(idom)%elems(ielem)%quad_pts(:)%c2_

        do igq = 1,size(x)


!   Two cylinder scattering
!            inA = ( x(igq) < -NINE + thickness )
!            inB = ( y(igq) >  FIVE - thickness )
!            inC = ( x(igq) >  NINE - thickness )
!            inD = ( y(igq) < -FIVE + thickness )
!
!
!            if ( inA ) then
!                fcn     =  abs( ( x - (-NINE+thickness) ) / thickness )**TWO
!                sigma_x = eps * fcn
!            
!            else if ( inC ) then
!                fcn     =  abs( ( x - (NINE-thickness) ) / thickness )**TWO
!                sigma_x = eps * fcn
!
!            else
!                sigma_x = ZERO
!
!            end if
!
!
!
!
!            if ( inB ) then
!                fcn     =  abs( ( y - ( FIVE-thickness) ) / thickness )**TWO
!                sigma_y = eps * fcn
!
!!            else if ( inD ) then
!!                fcn     =  abs( ( y - (-FIVE+thickness) ) / thickness )**TWO
!!                sigma_y = eps * fcn
!
!            else
!                sigma_y = ZERO
!            end if





            ! Munt duct
            inA = ( x(igq) < -THREE + thickness ) .and. ( y(igq) > 1.2_rk )
            inB = ( y(igq) >  4.6_rk - thickness )
            inC = ( x(igq) >  6.2_rk - thickness )
            inD = ( y(igq) < -FIVE + thickness )




            inA = .false.
            inB = .false.
            inC = .false.
            inD = .false.





            if ( inA ) then
                fcn(igq)     =  abs( ( x(igq) - (-THREE+thickness) ) / thickness )**TWO
                sigma_x(igq) = eps * fcn(igq)
            
            else if ( inC ) then
                fcn(igq)     =  abs( ( x(igq) - (6.2_rk-thickness) ) / thickness )**TWO
                sigma_x(igq) = eps * fcn(igq)

            else
                sigma_x(igq) = ZERO

            end if




            if ( inB ) then
                fcn(igq)     =  abs( ( y(igq) - ( 4.6_rk-thickness) ) / thickness )**TWO
                sigma_y(igq) = eps * fcn(igq)

!            else if ( inD ) then
!                fcn     =  abs( ( y - (-FIVE+thickness) ) / thickness )**TWO
!                sigma_y = eps * fcn

            else
                sigma_y(igq) = ZERO
            end if



            sigma(igq) = sigma_x(igq) * sigma_y(igq)


        end do









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


        !===========================
        !        MASS FLUX
        !===========================
        ! Unsteady source
        flux =            omega * rho_r

        ! F_r divergence source
!        flux = flux   +   (rhobar/y) * v_i

        ! F_theta divergence source
!        flux = flux   +   (rhobar*real(mod_m,rk)/y) * w_r


        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,irho_i,iblk,flux)


        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        ! Unsteady source
        flux =            omega * u_r
        
        ! F_r divergence source

        ! F_theta divergence source
        
        


        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,iu_i,iblk,flux)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        ! Unsteady source
        flux =            omega * v_r
        

        ! F_r divergence source
!        flux = flux     +   (ONE/(rhobar*y)) * p_i
        

        ! F_theta divergence source



        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,iv_i,iblk,flux)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        ! Unsteady source
        flux =            omega * w_r
        
        ! F_r divergence source


        ! F_theta divergence source
!        flux = flux   +   (real(mod_m,rk)/(rhobar*y)) * p_r
        

        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,iw_i,iblk,flux)

        !============================
        !       ENERGY FLUX
        !============================
        ! Unsteady source
        flux =            omega * p_r

        ! F_r divergence source
!        flux = flux   +   (gam*pbar/y)*v_i

        ! F_theta divergence source
!        flux = flux   +   (gam*pbar*real(mod_m,rk)/y)*w_r
        


        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,ip_i,iblk,flux)

    end subroutine compute
    !*********************************************************************************************************






end module PRIMLINEULERAXI_volume_advective_source_imag
