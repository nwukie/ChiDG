module PRIMLINEULER_volume_advective_sourceterms_real
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,FOUR,HALF,ZERO, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use type_mesh,              only: mesh_t
    use atype_volume_flux,      only: volume_flux_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    
    use mod_interpolate,        only: interpolate_element
    use mod_integrate,          only: integrate_volume_source
    use mod_DNAD_tools
    use DNAD_D

    use PRIMLINEULER_properties,    only: PRIMLINEULER_properties_t
    implicit none

    private


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    type, extends(volume_flux_t), public :: PRIMLINEULER_volume_advective_sourceterms_real_t


    contains

        procedure  :: compute
        
    end type PRIMLINEULER_volume_advective_sourceterms_real_t
    !*****************************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,idom,ielem,iblk)
        class(PRIMLINEULER_volume_advective_sourceterms_real_t),    intent(in)      :: self
        type(mesh_t),                                           intent(in)      :: mesh(:)
        type(solverdata_t),                                     intent(inout)   :: sdata
        class(properties_t),                                    intent(inout)   :: prop
        integer(ik),                                            intent(in)      :: idom, ielem, iblk

        ! Equation indices
        integer(ik)    :: irho_r,  irho_i
        integer(ik)    :: iu_r,    iu_i
        integer(ik)    :: iv_r,    iv_i
        integer(ik)    :: iw_r,    iw_i
        integer(ik)    :: ip_r,    ip_i


        integer(ik)    :: iseed, i, idonor, igq
        type(seed_t)   :: seed

        real(rk)    :: gam, omega, alpha, eps
        real(rk)    :: x, y, x0, y0



        type(AD_D), dimension(mesh(idom)%elems(ielem)%gq%vol%nnodes)      ::  &
                    rho_r, u_r, v_r, w_r, p_r, H,                        &
                    flux



        idonor = 0


        !
        ! Get equation indices
        !
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



        ! Initialize flux derivative storage
        flux = rho_r
        flux = ZERO

!        !===========================
!        !        MASS FLUX
!        !===========================
!        eps = ZERO
!        do igq = 1,size(rho)
!            x = mesh(idom)%elems(ielem)%quad_pts(igq)%c1_
!            y = mesh(idom)%elems(ielem)%quad_pts(igq)%c2_
!
!            flux(igq) = eps * exp(-alpha * (x**TWO + y**TWO) )
!
!        end do
!        flux = ZERO
!
!        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,irho_r,iblk,flux)
!
!
!        !===========================
!        !     X-MOMENTUM FLUX
!        !===========================
!        eps = ZERO
!        do igq = 1,size(rho)
!            x = mesh(idom)%elems(ielem)%quad_pts(igq)%c1_
!            y = mesh(idom)%elems(ielem)%quad_pts(igq)%c2_
!
!            flux(igq) = eps * exp(-alpha * (x**TWO + y**TWO) )
!
!        end do
!        flux = ZERO
!
!        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,irhou_r,iblk,flux)
!
!
!        !============================
!        !     Y-MOMENTUM FLUX
!        !============================
!        eps = ZERO
!        do igq = 1,size(rho)
!            x = mesh(idom)%elems(ielem)%quad_pts(igq)%c1_
!            y = mesh(idom)%elems(ielem)%quad_pts(igq)%c2_
!
!            flux(igq) = eps * exp(-alpha * (x**TWO + y**TWO) )
!
!        end do
!        flux = ZERO
!
!        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,irhov_r,iblk,flux)
!
!!        !============================
!!        !     Z-MOMENTUM FLUX
!!        !============================
!!        flux_x = (rhow*rhou)/rho
!!        flux_y = (rhow*rhov)/rho
!!        flux_z = (rhow*rhow)/rho  +  p
!!
!!        call integrate_volume_flux(mesh(idom)%elems(ielem),sdata,idom,irhow,iblk,flux_x,flux_y,flux_z)
!!
        !============================
        !       ENERGY FLUX
        !============================
        x0 = FOUR
        y0 = ZERO

        do igq = 1,size(rho_r)
            x = mesh(idom)%elems(ielem)%quad_pts(igq)%c1_
            y = mesh(idom)%elems(ielem)%quad_pts(igq)%c2_

            flux(igq) = exp(-LOG(TWO) * ((x-x0)**TWO + (y-y0)**TWO)/(0.2_rk**TWO) )

        end do

        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,ip_r,iblk,flux)

    end subroutine compute
    !**********************************************************************************************************






end module PRIMLINEULER_volume_advective_sourceterms_real
