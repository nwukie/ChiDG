module PRIMLINEULER_volume_advective_source_imag
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF,ZERO, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG,PI

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
    !-----------------------------------------------------------------------------------
    type, extends(volume_flux_t), public :: PRIMLINEULER_volume_advective_source_imag_t


    contains

        procedure  :: compute

    end type PRIMLINEULER_volume_advective_source_imag_t
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
        class(PRIMLINEULER_volume_advective_source_imag_t), intent(in)      :: self
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

        real(rk)    :: gam, omega, sigma_max, thickness, xl



        type(AD_D), dimension(mesh(idom)%elems(ielem)%gq%vol%nnodes)      ::    &
                    rho_r, u_r, v_r, w_r, p_r,                      & 
                    rho_i, u_i, v_i, w_i, p_i,                      &
                    p, H,                                                       &
                    flux

        real(rk), dimension(mesh(idom)%elems(ielem)%gq%vol%nnodes)      ::  &
                    x, sigma

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
        gam = 1.4_rk
        !omega = 348.329_rk * TWO * PI
        omega = 956._rk * TWO * PI



        !
        ! Compute PML coefficient
        !
        
        ! Get x-coordinate
        x = mesh(idom)%elems(ielem)%quad_pts(:)%c1_
        sigma_max = 600._rk
        thickness = 2._rk
        xl        = 8._rk

        sigma = sigma_max * abs( (x - xl)/thickness )**TWO


        do igq = 1,size(x)
            if ( x(igq) < xl ) then
                sigma(igq) = ZERO
            end if
        end do

        sigma = ZERO

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
        flux =  omega * rho_r  -  sigma*rho_i

        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,irho_i,iblk,flux)


        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux =  omega * u_r  -  sigma*u_i

        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,iu_i,iblk,flux)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux =  omega * v_r  -  sigma*v_i

        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,iv_i,iblk,flux)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        flux =  omega * w_r  -  sigma*w_i

        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,iw_i,iblk,flux)

        !============================
        !       ENERGY FLUX
        !============================
        flux =  omega * p_r  -  sigma*p_i

        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,ip_i,iblk,flux)

    end subroutine compute
    !*********************************************************************************************************






end module PRIMLINEULER_volume_advective_source_imag
