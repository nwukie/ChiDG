module PRIMLINEULER_LaxFriedrichs_flux_imag
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ONE,TWO,HALF,FOUR,ZERO, ME, NEIGHBOR

    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D

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
    type, extends(operator_t), public :: PRIMLINEULER_LaxFriedrichs_flux_imag_t

    contains

        procedure   :: init
        procedure   :: compute

    end type PRIMLINEULER_LaxFriedrichs_flux_imag_t
    !********************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2018
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init(self)
        class(PRIMLINEULER_LaxFriedrichs_flux_imag_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('PRIMLINEULER LaxFriedrichs Imag')

        ! Set operator type
        call self%set_operator_type('Boundary Advective Operator')

        ! Set operator equations
        call self%add_primary_field("Density(imag)"   )
        call self%add_primary_field("Velocity-1(imag)")
        call self%add_primary_field("Velocity-2(imag)")
        call self%add_primary_field("Velocity-3(imag)")
        call self%add_primary_field("Pressure(imag)"  )


    end subroutine init
    !*****************************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2018
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
            un, wave_c, unorm_1, unorm_2, unorm_3


        !
        ! Interpolate solution to quadrature nodes
        !
        rho_m = worker%get_field('Density(imag)',    'value', 'face interior')
        rho_p = worker%get_field('Density(imag)',    'value', 'face exterior')

        u_m   = worker%get_field('Velocity-1(imag)', 'value', 'face interior')
        u_p   = worker%get_field('Velocity-1(imag)', 'value', 'face exterior')

        v_m   = worker%get_field('Velocity-2(imag)', 'value', 'face interior')
        v_p   = worker%get_field('Velocity-2(imag)', 'value', 'face exterior')

        w_m   = worker%get_field('Velocity-3(imag)', 'value', 'face interior')
        w_p   = worker%get_field('Velocity-3(imag)', 'value', 'face exterior')

        p_m   = worker%get_field('Pressure(imag)',   'value', 'face interior')
        p_p   = worker%get_field('Pressure(imag)',   'value', 'face exterior')


        !
        ! Get normal vector
        !
        unorm_1 = worker%unit_normal_ale(1)
        unorm_2 = worker%unit_normal_ale(2)
        unorm_3 = worker%unit_normal_ale(3)


!        wave = rho_m
!        do igq = 1,size(wave)
!            wave(igq)%x_ad_  = ZERO
!            wave(igq)%xp_ad_ = ZERO
!        end do



        !--------------------------------------
        !  Compute wave speeds
        !--------------------------------------



        !
        ! Compute normal velocities: dot-product vector projection along unit-normal direction
        !
        !un = unormx*ubar + unormy*vbar + unormz*wbar
        un = unorm_1*ubar + unorm_2*vbar + unorm_3*wbar

        !
        ! Compute wave speeds
        !
        !wave_c = abs(un) + cbar
        !wave = wave_c
        wave= abs(un) + cbar





        !================================
        !       MASS FLUX
        !================================
        !upwind = -wave*(rho_p - rho_m)
        !integrand = HALF*(upwind*normx*unormx  +  upwind*normy*unormy  +  upwind*normz*unormz )

        dissipation = HALF*wave*(rho_m - rho_p)

        call worker%integrate_boundary_upwind('Density(imag)',dissipation)


        !================================
        !       X-MOMENTUM FLUX
        !================================
        !upwind = -wave*(u_p - u_m)
        !integrand = HALF*(upwind*normx*unormx  +  upwind*normy*unormy  +  upwind*normz*unormz )
        !call worker%integrate_boundary(iu,integrand)

        dissipation = HALF*wave*(u_m - u_p)

        call worker%integrate_boundary_upwind('Velocity-1(imag)',dissipation)


        !================================
        !       Y-MOMENTUM FLUX
        !================================
        !upwind = -wave*(v_p - v_m)
        !integrand = HALF*(upwind*normx*unormx  +  upwind*normy*unormy  +  upwind*normz*unormz )
        !call worker%integrate_boundary(iv,integrand)

        dissipation = HALF*wave*(v_m - v_p)

        call worker%integrate_boundary_upwind('Velocity-2(imag)',dissipation)

        !================================
        !       Z-MOMENTUM FLUX
        !================================
        !upwind = -wave*(w_p - w_m)
        !integrand = HALF*(upwind*normx*unormx  +  upwind*normy*unormy  +  upwind*normz*unormz )
        !call worker%integrate_boundary(iw,integrand)

        dissipation = HALF*wave*(w_m - w_p)

        call worker%integrate_boundary_upwind('Velocity-3(imag)',dissipation)


        !================================
        !          ENERGY FLUX
        !================================
        !upwind = -wave*(p_p - p_m)
        !integrand = HALF*(upwind*normx*unormx  +  upwind*normy*unormy  +  upwind*normz*unormz )
        !call worker%integrate_boundary(ip,integrand)

        dissipation = HALF*wave*(p_m - p_p)

        call worker%integrate_boundary_upwind('Pressure(imag)',dissipation)



    end subroutine compute
    !*****************************************************************************************













end module PRIMLINEULER_LaxFriedrichs_flux_imag
