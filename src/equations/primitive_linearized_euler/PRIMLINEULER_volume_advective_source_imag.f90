module PRIMLINEULER_volume_advective_source_imag
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,THREE,FOUR,FIVE,EIGHT,NINE,HALF,ZERO

    use type_volume_flux,       only: volume_flux_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
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
    subroutine compute(self,worker,prop)
        class(PRIMLINEULER_volume_advective_source_imag_t), intent(in)      :: self
        type(chidg_worker_t),                               intent(inout)   :: worker
        class(properties_t),                                intent(inout)   :: prop

        ! Equation indices
        integer(ik)    :: irho_r,  irho_i
        integer(ik)    :: iu_r, iu_i
        integer(ik)    :: iv_r, iv_i
        integer(ik)    :: iw_r, iw_i
        integer(ik)    :: ip_r, ip_i

        integer(ik) :: igq


        type(AD_D), allocatable, dimension(:)   ::  &
            rho_r, u_r, v_r, w_r, p_r,              & 
            rho_i, u_i, v_i, w_i, p_i,              &
            p, H,                                   &
            source

        real(rk),   allocatable, dimension(:)   ::  &
            x, y, z, r, sigma_x, sigma_y, sigma_z, fcn

        logical :: inA = .false.
        logical :: inB = .false.
        logical :: inC = .false.
        logical :: inD = .false.
        logical :: inE = .false.
        logical :: inF = .false.



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



        !===========================
        !        MASS FLUX
        !===========================
        source =  omega * rho_r * (ONE - (sigma_x*sigma_y+sigma_x*sigma_z+sigma_y*sigma_z)/(omega*omega) )  + &
                  omega * rho_i * ( -(sigma_x+sigma_y+sigma_z)/omega  +  (sigma_x*sigma_y*sigma_z)/(omega*omega*omega) )

        call worker%integrate_volume(irho_i,source)


        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        source =  omega * u_r * (ONE - (sigma_x*sigma_y+sigma_x*sigma_z+sigma_y*sigma_z)/(omega*omega) )    + &
                  omega * u_i * ( -(sigma_x+sigma_y+sigma_z)/omega  +  (sigma_x*sigma_y*sigma_z)/(omega*omega*omega) )

        call worker%integrate_volume(iu_i,source)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        source =  omega * v_r * (ONE - (sigma_x*sigma_y+sigma_x*sigma_z+sigma_y*sigma_z)/(omega*omega) )    + &
                  omega * v_i * ( -(sigma_x+sigma_y+sigma_z)/omega  +  (sigma_x*sigma_y*sigma_z)/(omega*omega*omega) )

        call worker%integrate_volume(iv_i,source)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        source =  omega * w_r * (ONE - (sigma_x*sigma_y+sigma_x*sigma_z+sigma_y*sigma_z)/(omega*omega) )    + &
                  omega * w_i * ( -(sigma_x+sigma_y+sigma_z)/omega  +  (sigma_x*sigma_y*sigma_z)/(omega*omega*omega) )

        call worker%integrate_volume(iw_i,source)

        !============================
        !       ENERGY FLUX
        !============================
        source =  omega * p_r * (ONE - (sigma_x*sigma_y+sigma_x*sigma_z+sigma_y*sigma_z)/(omega*omega) )    + &
                  omega * p_i * ( -(sigma_x+sigma_y+sigma_z)/omega  +  (sigma_x*sigma_y*sigma_z)/(omega*omega*omega) )

        call worker%integrate_volume(ip_i,source)

    end subroutine compute
    !*********************************************************************************************************






end module PRIMLINEULER_volume_advective_source_imag
