module PRIMLINEULER_LaxFriedrichs_flux_imag
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF,FOUR,ZERO, ME, NEIGHBOR

    use type_boundary_flux,     only: boundary_flux_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D

    use PRIMLINEULER_properties,    only: PRIMLINEULER_properties_t
    use mod_primitive_linearized_euler
    implicit none

    private


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/17/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------
    type, extends(boundary_flux_t), public :: PRIMLINEULER_LaxFriedrichs_flux_imag_t

    contains

        procedure  :: compute

    end type PRIMLINEULER_LaxFriedrichs_flux_imag_t
    !********************************************************************************************










contains






    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/16/2017
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(PRIMLINEULER_LaxFriedrichs_flux_imag_t),  intent(in)      :: self
        type(chidg_worker_t),                           intent(inout)   :: worker
        class(properties_t),                            intent(inout)   :: prop

        ! Equation indices
        integer(ik) :: irho
        integer(ik) :: iu
        integer(ik) :: iv
        integer(ik) :: iw
        integer(ik) :: ip

        integer(ik) :: igq


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            rho_m,   rho_p,                         &
            u_m,     u_p,                           &
            v_m,     v_p,                           &
            w_m,     w_p,                           &
            p_m,     p_p,                           &
            integrand,  upwind,                     &
            wave 

        real(rk),   allocatable, dimension(:)   ::  &
            un, wave_c, normx, normy, normz, unormx, unormy, unormz


        irho = prop%get_eqn_index("rho_i")
        iu   = prop%get_eqn_index("u_i")
        iv   = prop%get_eqn_index("v_i")
        iw   = prop%get_eqn_index("w_i")
        ip   = prop%get_eqn_index("p_i")



        !
        ! Interpolate solution to quadrature nodes
        !
        rho_m  = worker%interpolate(irho, 'value', ME)
        rho_p  = worker%interpolate(irho, 'value', NEIGHBOR)

        u_m    = worker%interpolate(iu,   'value', ME)
        u_p    = worker%interpolate(iu,   'value', NEIGHBOR)

        v_m    = worker%interpolate(iv,   'value', ME)
        v_p    = worker%interpolate(iv,   'value', NEIGHBOR)

        w_m    = worker%interpolate(iw,   'value', ME)
        w_p    = worker%interpolate(iw,   'value', NEIGHBOR)

        p_m    = worker%interpolate(ip,   'value', ME)
        p_p    = worker%interpolate(ip,   'value', NEIGHBOR)




        normx  = worker%normal(1)
        normy  = worker%normal(2)
        normz  = worker%normal(3)

        unormx = worker%unit_normal(1)
        unormy = worker%unit_normal(2)
        unormz = worker%unit_normal(3)


        wave = rho_m
        do igq = 1,size(wave)
            wave(igq)%x_ad_  = ZERO
            wave(igq)%xp_ad_ = ZERO
        end do



        !--------------------------------------
        !  Compute wave speeds
        !--------------------------------------



        !
        ! Compute normal velocities: dot-product vector projection along unit-normal direction
        !
        un = unormx*ubar + unormy*vbar + unormz*wbar

        !
        ! Compute wave speeds
        !
        wave_c = abs(un) + cbar

        wave = wave_c




        !================================
        !       MASS FLUX
        !================================
        upwind = -wave*(rho_p - rho_m)

        integrand = HALF*(upwind*normx*unormx  +  upwind*normy*unormy  +  upwind*normz*unormz )

        call worker%integrate_boundary(irho,integrand)


        !================================
        !       X-MOMENTUM FLUX
        !================================
        upwind = -wave*(u_p - u_m)

        integrand = HALF*(upwind*normx*unormx  +  upwind*normy*unormy  +  upwind*normz*unormz )

        call worker%integrate_boundary(iu,integrand)


        !================================
        !       Y-MOMENTUM FLUX
        !================================
        upwind = -wave*(v_p - v_m)

        integrand = HALF*(upwind*normx*unormx  +  upwind*normy*unormy  +  upwind*normz*unormz )

        call worker%integrate_boundary(iv,integrand)

        !================================
        !       Z-MOMENTUM FLUX
        !================================
        upwind = -wave*(w_p - w_m)

        integrand = HALF*(upwind*normx*unormx  +  upwind*normy*unormy  +  upwind*normz*unormz )

        call worker%integrate_boundary(iw,integrand)

        !================================
        !          ENERGY FLUX
        !================================
        upwind = -wave*(p_p - p_m)

        integrand = HALF*(upwind*normx*unormx  +  upwind*normy*unormy  +  upwind*normz*unormz )

        call worker%integrate_boundary(ip,integrand)


    end subroutine compute
    !**********************************************************************************************************************













end module PRIMLINEULER_LaxFriedrichs_flux_imag
