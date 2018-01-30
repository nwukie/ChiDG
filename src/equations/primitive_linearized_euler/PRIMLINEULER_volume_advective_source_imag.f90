module PRIMLINEULER_volume_advective_source_imag
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ONE,TWO,THREE,FOUR,FIVE,EIGHT,NINE,HALF,ZERO

    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
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
    type, extends(operator_t), public :: PRIMLINEULER_volume_advective_source_imag_t


    contains

        procedure   :: init
        procedure   :: compute

    end type PRIMLINEULER_volume_advective_source_imag_t
    !************************************************************************************










contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(PRIMLINEULER_volume_advective_source_imag_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("PRIMLINEULER Volume Source Imag")

        ! Set operator type
        call self%set_operator_type("Volume Advective Flux")

        ! Set operator equations
        call self%add_primary_field("Density(imag)"   )
        call self%add_primary_field("Velocity-1(imag)")
        call self%add_primary_field("Velocity-2(imag)")
        call self%add_primary_field("Velocity-3(imag)")
        call self%add_primary_field("Pressure(imag)"  )

    end subroutine init
    !********************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2018
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



!        !
!        ! Get coordinates
!        !
!        x = worker%x('boundary')
!        y = worker%y('boundary')
!        z = worker%z('boundary')
!        r = sqrt(y**TWO + z**TWO)
!
!        !
!        ! Compute PML Layers
!        !
!        do igq = 1,size(x)
!
!
!            ! Munt duct
!            inA = ( x(igq) < -THREE + thickness ) .and. ( r(igq) > 1.212_rk )
!            inB = ( x(igq) >  THREE - thickness )
!            inC = ( y(igq) < -2.121_rk + thickness )
!            inD = ( y(igq) >  2.121_rk - thickness )
!            inE = ( z(igq) < -2.121_rk + thickness )
!            inF = ( z(igq) >  2.121_rk - thickness )
!
!!            ! Monopole
!!            inA = ( x(igq) < -100._rk + thickness ) 
!!            inB = ( x(igq) >  100._rk - thickness )
!!            inC = ( y(igq) < -100._rk + thickness )
!!            inD = ( y(igq) >  100._rk - thickness )
!!            inE = ( z(igq) < -2.121_rk + thickness )
!!            inF = ( z(igq) >  2.121_rk - thickness )
!
!
!!            inA = .false.
!!            inB = .false.
!!            inC = .false.
!!            inD = .false.
!!            inE = .false.
!!            inF = .false.
!
!
!            ! X-PML
!            if ( inA ) then
!                fcn(igq)     =  abs( ( x(igq) - (-3._rk+thickness) ) / thickness )**TWO
!                sigma_x(igq) = eps * fcn(igq)
!            else if ( inB ) then
!                fcn(igq)     =  abs( ( x(igq) - ( 3._rk-thickness) ) / thickness )**TWO
!                sigma_x(igq) = eps * fcn(igq)
!            else
!                sigma_x(igq) = ZERO
!            end if
!
!
!            ! Y-PML
!            if ( inC ) then
!                fcn(igq)     =  abs( ( y(igq) - (-2.121_rk+thickness) ) / thickness )**TWO
!                sigma_y(igq) = eps * fcn(igq)
!            else if ( inD ) then
!                fcn(igq)     =  abs( ( y(igq) - ( 2.121_rk-thickness) ) / thickness )**TWO
!                sigma_y(igq) = eps * fcn(igq)
!            else
!                sigma_y(igq) = ZERO
!            end if
!
!
!            ! Z-PML
!            if ( inE ) then
!                fcn(igq)     =  abs( ( z(igq) - (-2.121_rk+thickness) ) / thickness )**TWO
!                sigma_z(igq) = eps * fcn(igq)
!            else if ( inF ) then
!                fcn(igq)     =  abs( ( z(igq) - ( 2.121_rk-thickness) ) / thickness )**TWO
!                sigma_z(igq) = eps * fcn(igq)
!            else
!                sigma_z(igq) = ZERO
!            end if
!
!        end do





        !
        ! Interpolate solution to quadrature nodes: REAL
        !
        rho_r = worker%get_field('Density(real)'   , 'value', 'element')
        u_r   = worker%get_field('Velocity-1(real)', 'value', 'element')
        v_r   = worker%get_field('Velocity-2(real)', 'value', 'element')
        w_r   = worker%get_field('Velocity-3(real)', 'value', 'element')
        p_r   = worker%get_field('Pressure(real)'  , 'value', 'element')

        !
        ! Interpolate solution to quadrature nodes: IMAG
        !
        rho_i = worker%get_field('Density(imag)'   , 'value', 'element')
        u_i   = worker%get_field('Velocity-1(imag)', 'value', 'element')
        v_i   = worker%get_field('Velocity-2(imag)', 'value', 'element')
        w_i   = worker%get_field('Velocity-3(imag)', 'value', 'element')
        p_i   = worker%get_field('Pressure(imag)'  , 'value', 'element')







        !===========================
        !        MASS FLUX
        !===========================
        source =  omega * rho_r 

        call worker%integrate_volume_source('Density(imag)',source)


        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        source =  omega * u_r 

        call worker%integrate_volume_source('Velocity-1(imag)',source)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        source =  omega * v_r 

        call worker%integrate_volume_source('Velocity-2(imag)',source)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        source =  omega * w_r 

        call worker%integrate_volume_source('Velocity-3(imag)',source)

        !============================
        !       ENERGY FLUX
        !============================
        source =  omega * p_r 

        call worker%integrate_volume_source('Pressure(imag)',source)

    end subroutine compute
    !*********************************************************************************************************






end module PRIMLINEULER_volume_advective_source_imag
