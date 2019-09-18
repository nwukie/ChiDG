module PRIMLINEULER_volume_advective_flux_imag
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ONE,TWO,THREE,FOUR,FIVE,EIGHT,NINE,HALF,ZERO, PI

    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D

    use mod_primitive_linearized_euler
    implicit none
    private


    !>  Volume advective flux for Linearized Euler equations - imaginary.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2018
    !!
    !!
    !----------------------------------------------------------------------------------
    type, extends(operator_t), public :: PRIMLINEULER_volume_advective_flux_imag_t


    contains

        procedure   :: init
        procedure   :: compute

    end type PRIMLINEULER_volume_advective_flux_imag_t
    !***********************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(PRIMLINEULER_volume_advective_flux_imag_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("PRIMLINEULER Volume Flux Imag")

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
    !----------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(PRIMLINEULER_volume_advective_flux_imag_t),   intent(in)      :: self
        type(chidg_worker_t),                               intent(inout)   :: worker
        class(properties_t),                                intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:)   ::  &
            rho_r, u_r, v_r, w_r, p_r,              &
            rho_i, u_i, v_i, w_i, p_i,              &
            flux_1, flux_2, flux_3

        real(rk),   allocatable, dimension(:)   ::  &
            x, y, z, r, sigma_x, sigma_y, sigma_z, fcn

        logical :: inA = .false.
        logical :: inB = .false.
        logical :: inC = .false.
        logical :: inD = .false.
        logical :: inE = .false.
        logical :: inF = .false.



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
!
!!            ! Two-cylinder scattering
!!            inA = ( x(igq) < -NINE + thickness )
!!            inB = ( x(igq) >  NINE - thickness )
!!            inC = ( y(igq) < -2.121_rk + thickness )
!!            inD = ( y(igq) >  FIVE - thickness )
!!            inE = ( z(igq) < -2.121_rk + thickness )
!!            inF = ( z(igq) >  2.121_rk - thickness )
!
!
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








        !===========================
        !        MASS FLUX
        !===========================
        flux_1 = rho_1_rho * rho_i  + &
                 rho_1_u   * u_i    + &
                 rho_1_v   * v_i    + &
                 rho_1_w   * w_i    + &
                 rho_1_p   * p_i

        flux_2 = rho_2_rho * rho_i  + &
                 rho_2_u   * u_i    + &
                 rho_2_v   * v_i    + &
                 rho_2_w   * w_i    + &
                 rho_2_p   * p_i

        flux_3 = rho_3_rho * rho_i  + &
                 rho_3_u   * u_i    + &
                 rho_3_v   * v_i    + &
                 rho_3_w   * w_i    + &
                 rho_3_p   * p_i

        call worker%integrate_volume_flux('Density(imag)','Advection',flux_1,flux_2,flux_3)



        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux_1 = u_1_rho * rho_i  + &
                 u_1_u   * u_i    + &
                 u_1_v   * v_i    + &
                 u_1_w   * w_i    + &
                 u_1_p   * p_i

        flux_2 = u_2_rho * rho_i  + &
                 u_2_u   * u_i    + &
                 u_2_v   * v_i    + &
                 u_2_w   * w_i    + &
                 u_2_p   * p_i

        flux_3 = u_3_rho * rho_i  + &
                 u_3_u   * u_i    + &
                 u_3_v   * v_i    + &
                 u_3_w   * w_i    + &
                 u_3_p   * p_i

        call worker%integrate_volume_flux('Velocity-1(imag)','Advection',flux_1,flux_2,flux_3)



        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux_1 = v_1_rho * rho_i  + &
                 v_1_u   * u_i    + &
                 v_1_v   * v_i    + &
                 v_1_w   * w_i    + &
                 v_1_p   * p_i

        flux_2 = v_2_rho * rho_i  + &
                 v_2_u   * u_i    + &
                 v_2_v   * v_i    + &
                 v_2_w   * w_i    + &
                 v_2_p   * p_i

        flux_3 = v_3_rho * rho_i  + &
                 v_3_u   * u_i    + &
                 v_3_v   * v_i    + &
                 v_3_w   * w_i    + &
                 v_3_p   * p_i

        call worker%integrate_volume_flux('Velocity-2(imag)','Advection',flux_1,flux_2,flux_3)


        !============================
        !     Z-MOMENTUM FLUX
        !============================
        flux_1 = w_1_rho * rho_i  + &
                 w_1_u   * u_i    + &
                 w_1_v   * v_i    + &
                 w_1_w   * w_i    + &
                 w_1_p   * p_i

        flux_2 = w_2_rho * rho_i  + &
                 w_2_u   * u_i    + &
                 w_2_v   * v_i    + &
                 w_2_w   * w_i    + &
                 w_2_p   * p_i

        flux_3 = w_3_rho * rho_i  + &
                 w_3_u   * u_i    + &
                 w_3_v   * v_i    + &
                 w_3_w   * w_i    + &
                 w_3_p   * p_i

        call worker%integrate_volume_flux('Velocity-3(imag)','Advection',flux_1,flux_2,flux_3)


        !============================
        !       ENERGY FLUX
        !============================
        flux_1 = p_1_rho * rho_i  + &
                 p_1_u   * u_i    + &
                 p_1_v   * v_i    + &
                 p_1_w   * w_i    + &
                 p_1_p   * p_i

        flux_2 = p_2_rho * rho_i  + &
                 p_2_u   * u_i    + &
                 p_2_v   * v_i    + &
                 p_2_w   * w_i    + &
                 p_2_p   * p_i
  
        flux_3 = p_3_rho * rho_i  + &
                 p_3_u   * u_i    + &
                 p_3_v   * v_i    + &
                 p_3_w   * w_i    + &
                 p_3_p   * p_i

        call worker%integrate_volume_flux('Pressure(imag)','Advection',flux_1,flux_2,flux_3)


    end subroutine compute
    !******************************************************************************************************






end module PRIMLINEULER_volume_advective_flux_imag
