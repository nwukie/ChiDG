module LINEULER_volume_advective_source_real
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ONE,TWO,HALF,ZERO, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG,PI

    use type_mesh,              only: mesh_t
    use type_volume_flux,       only: volume_flux_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    use type_element_info,      only: element_info_t
    use type_function_info,     only: function_info_t
    
    use mod_interpolate,        only: interpolate_element
    use mod_integrate,          only: integrate_volume_source
    use DNAD_D

    use LINEULER_properties,    only: LINEULER_properties_t
    implicit none

    private


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    type, extends(volume_flux_t), public :: LINEULER_volume_advective_source_real_t


    contains

        procedure  :: compute

    end type LINEULER_volume_advective_source_real_t
    !***************************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine compute(self,mesh,sdata,prop,elem_info,function_info)
        class(LINEULER_volume_advective_source_real_t), intent(in)      :: self
        type(mesh_t),                                   intent(in)      :: mesh(:)
        type(solverdata_t),                             intent(inout)   :: sdata
        class(properties_t),                            intent(inout)   :: prop
        type(element_info_t),                           intent(in)      :: elem_info
        type(function_info_t),                          intent(in)      :: function_info

        ! Equation indices
        integer(ik)    :: irho_r, irho_i
        integer(ik)    :: irhou_r, irhou_i
        integer(ik)    :: irhov_r, irhov_i
        integer(ik)    :: irhow_r, irhow_i
        integer(ik)    :: irhoE_r, irhoE_i


        real(rk)       :: gam, omega, sigma_max, thickness, xl
        integer(ik)    :: idom, ielem, iblk, igq



        type(AD_D), dimension(mesh(elem_info%idomain_l)%elems(elem_info%ielement_l)%gq%vol%nnodes)      ::    &
                    rho_r, rhou_r, rhov_r, rhow_r, rhoE_r,                      &
                    rho_i, rhou_i, rhov_i, rhow_i, rhoE_i,                      &
                    p,     H,                                                   &
                    flux

        real(rk), dimension(mesh(elem_info%idomain_l)%elems(elem_info%ielement_l)%gq%vol%nnodes)      ::  &
                    x, sigma

        idom  = elem_info%idomain_l
        ielem = elem_info%ielement_l
        iblk  = function_info%iblk


        !-------------------------------------------------------------
        irho_r  = prop%get_eqn_index("rho_r")
        irhou_r = prop%get_eqn_index("rhou_r")
        irhov_r = prop%get_eqn_index("rhov_r")
        irhow_r = prop%get_eqn_index("rhow_r")
        irhoE_r = prop%get_eqn_index("rhoE_r")

        irho_i  = prop%get_eqn_index("rho_i")
        irhou_i = prop%get_eqn_index("rhou_i")
        irhov_i = prop%get_eqn_index("rhov_i")
        irhow_i = prop%get_eqn_index("rhow_i")
        irhoE_i = prop%get_eqn_index("rhoE_i")



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
        ! Interpolate solution to quadrature nodes
        !
        call interpolate_element(mesh,sdata%q,idom,ielem,irho_i, rho_i, function_info%seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhou_i,rhou_i,function_info%seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhov_i,rhov_i,function_info%seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhow_i,rhow_i,function_info%seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhoE_i,rhoE_i,function_info%seed)

        call interpolate_element(mesh,sdata%q,idom,ielem,irho_r, rho_r, function_info%seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhou_r,rhou_r,function_info%seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhov_r,rhov_r,function_info%seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhow_r,rhow_r,function_info%seed)
        call interpolate_element(mesh,sdata%q,idom,ielem,irhoE_r,rhoE_r,function_info%seed)


        !===========================
        !        MASS FLUX
        !===========================
        flux = -omega * rho_i   -  sigma*rho_r

        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,irho_r,iblk,flux)


        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux = -omega * rhou_i   -  sigma*rhou_r

        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,irhou_r,iblk,flux)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux = -omega * rhov_i   -  sigma*rhov_r

        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,irhov_r,iblk,flux)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        flux = -omega * rhow_i   -  sigma*rhow_r

        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,irhow_r,iblk,flux)

        !============================
        !       ENERGY FLUX
        !============================
        flux = -omega * rhoE_i   -  sigma*rhoE_r

        call integrate_volume_source(mesh(idom)%elems(ielem),sdata,idom,irhoE_r,iblk,flux)

    end subroutine compute
    !*********************************************************************************************************






end module LINEULER_volume_advective_source_real
