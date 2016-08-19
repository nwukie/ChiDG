module PRIMLINEULER_volume_advective_source_imag
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,THREE,FOUR,FIVE,EIGHT,NINE,HALF,ZERO

    use type_mesh,              only: mesh_t
    use type_volume_flux,       only: volume_flux_t
    use type_solverdata,        only: solverdata_t
    use type_properties,        only: properties_t
    use type_element_info,      only: element_info_t
    use type_function_info,     only: function_info_t
    
    use mod_interpolate,        only: interpolate
    use mod_integrate,          only: integrate_volume_source
    use DNAD_D

    use PRIMLINEULER_properties,        only: PRIMLINEULER_properties_t
    use mod_primitive_linearized_euler, only: omega, gam, thickness, eps, rhobar, pbar, mod_m
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
    subroutine compute(self,mesh,sdata,prop,elem_info,function_info)
        class(PRIMLINEULER_volume_advective_source_imag_t), intent(in)      :: self
        type(mesh_t),                                   intent(in)      :: mesh(:)
        type(solverdata_t),                             intent(inout)   :: sdata
        class(properties_t),                            intent(inout)   :: prop
        type(element_info_t),                           intent(in)      :: elem_info
        type(function_info_t),                          intent(in)      :: function_info

        ! Equation indices
        integer(ik)    :: irho_r,  irho_i
        integer(ik)    :: iu_r, iu_i
        integer(ik)    :: iv_r, iv_i
        integer(ik)    :: iw_r, iw_i
        integer(ik)    :: ip_r, ip_i

        integer(ik)    :: idom, ielem, igq


        type(AD_D), dimension(mesh(elem_info%idomain_l)%elems(elem_info%ielement_l)%gq%vol%nnodes)      ::  &
                    rho_r, u_r, v_r, w_r, p_r,                      & 
                    rho_i, u_i, v_i, w_i, p_i,                      &
                    p, H,                                           &
                    flux

        real(rk), dimension(mesh(elem_info%idomain_l)%elems(elem_info%ielement_l)%gq%vol%nnodes)      ::  &
                    x, y, z, r, sigma_x, sigma_y, sigma_z, fcn

        logical :: inA = .false.
        logical :: inB = .false.
        logical :: inC = .false.
        logical :: inD = .false.
        logical :: inE = .false.
        logical :: inF = .false.

        idom  = elem_info%idomain_l
        ielem = elem_info%ielement_l


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



        !
        ! Get coordinates
        !
        x = mesh(idom)%elems(ielem)%quad_pts(:)%c1_
        y = mesh(idom)%elems(ielem)%quad_pts(:)%c2_
        z = mesh(idom)%elems(ielem)%quad_pts(:)%c3_
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
!            inA = ( x(igq) < -100._rk + thickness ) 
!            inB = ( x(igq) >  100._rk - thickness )
!            inC = ( y(igq) < -100._rk + thickness )
!            inD = ( y(igq) >  100._rk - thickness )
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



        !
        ! Interpolate solution to quadrature nodes
        !
        rho_r = interpolate(mesh,sdata,elem_info,function_info,irho_r, 'value')
        u_r   = interpolate(mesh,sdata,elem_info,function_info,iu_r,   'value')
        v_r   = interpolate(mesh,sdata,elem_info,function_info,iv_r,   'value')
        w_r   = interpolate(mesh,sdata,elem_info,function_info,iw_r,   'value')
        p_r   = interpolate(mesh,sdata,elem_info,function_info,ip_r,   'value')

        rho_i = interpolate(mesh,sdata,elem_info,function_info,irho_i, 'value')
        u_i   = interpolate(mesh,sdata,elem_info,function_info,iu_i,   'value')
        v_i   = interpolate(mesh,sdata,elem_info,function_info,iv_i,   'value')
        w_i   = interpolate(mesh,sdata,elem_info,function_info,iw_i,   'value')
        p_i   = interpolate(mesh,sdata,elem_info,function_info,ip_i,   'value')

        !===========================
        !        MASS FLUX
        !===========================
        flux =  omega * rho_r * (ONE - (sigma_x*sigma_y+sigma_x*sigma_z+sigma_y*sigma_z)/(omega*omega) )  + &
                omega * rho_i * ( -(sigma_x+sigma_y+sigma_z)/omega  +  (sigma_x*sigma_y*sigma_z)/(omega*omega*omega) )

        call integrate_volume_source(mesh,sdata,elem_info,function_info,irho_i,flux)


        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux =  omega * u_r * (ONE - (sigma_x*sigma_y+sigma_x*sigma_z+sigma_y*sigma_z)/(omega*omega) )    + &
                omega * u_i * ( -(sigma_x+sigma_y+sigma_z)/omega  +  (sigma_x*sigma_y*sigma_z)/(omega*omega*omega) )

        call integrate_volume_source(mesh,sdata,elem_info,function_info,iu_i,flux)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux =  omega * v_r * (ONE - (sigma_x*sigma_y+sigma_x*sigma_z+sigma_y*sigma_z)/(omega*omega) )    + &
                omega * v_i * ( -(sigma_x+sigma_y+sigma_z)/omega  +  (sigma_x*sigma_y*sigma_z)/(omega*omega*omega) )

        call integrate_volume_source(mesh,sdata,elem_info,function_info,iv_i,flux)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        flux =  omega * w_r * (ONE - (sigma_x*sigma_y+sigma_x*sigma_z+sigma_y*sigma_z)/(omega*omega) )    + &
                omega * w_i * ( -(sigma_x+sigma_y+sigma_z)/omega  +  (sigma_x*sigma_y*sigma_z)/(omega*omega*omega) )

        call integrate_volume_source(mesh,sdata,elem_info,function_info,iw_i,flux)

        !============================
        !       ENERGY FLUX
        !============================
        flux =  omega * p_r * (ONE - (sigma_x*sigma_y+sigma_x*sigma_z+sigma_y*sigma_z)/(omega*omega) )    + &
                omega * p_i * ( -(sigma_x+sigma_y+sigma_z)/omega  +  (sigma_x*sigma_y*sigma_z)/(omega*omega*omega) )

        call integrate_volume_source(mesh,sdata,elem_info,function_info,ip_i,flux)

    end subroutine compute
    !*********************************************************************************************************






end module PRIMLINEULER_volume_advective_source_imag
