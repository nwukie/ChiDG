module PRIMLINEULER_LaxFriedrichs_flux_real
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF,FOUR,ZERO,ME, NEIGHBOR

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
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
    !-------------------------------------------------------------------------------------
    type, extends(operator_t), public :: PRIMLINEULER_LaxFriedrichs_flux_real_t

    contains

        procedure   :: init
        procedure   :: compute

    end type PRIMLINEULER_LaxFriedrichs_flux_real_t
    !**************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2018
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init(self)
        class(PRIMLINEULER_LaxFriedrichs_flux_real_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('PRIMLINEULER LaxFriedrichs Real')

        ! Set operator type
        call self%set_operator_type('Boundary Advective Operator')

        ! Set operator equations
        call self%add_primary_field("Density(real)"   )
        call self%add_primary_field("Velocity-1(real)")
        call self%add_primary_field("Velocity-2(real)")
        call self%add_primary_field("Velocity-3(real)")
        call self%add_primary_field("Pressure(real)"  )


    end subroutine init
    !*****************************************************************************************
    


    !>  Real component of numerical flux dissipation
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2018
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(PRIMLINEULER_LaxFriedrichs_flux_real_t),  intent(in)      :: self
        type(chidg_worker_t),                           intent(inout)   :: worker
        class(properties_t),                            intent(inout)   :: prop



        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            rho_m,   rho_p,                         &
            u_m,     u_p,                           &
            v_m,     v_p,                           &
            w_m,     w_p,                           &
            p_m,     p_p,                           &
            dissipation,  wave


        real(rk),   allocatable, dimension(:)   ::  &
            un, wave_c, unorm_1, unorm_2, unorm_3


        !
        ! Interpolate solution to quadrature nodes
        !
        rho_m = worker%get_field('Density(real)',    'value', 'face interior')
        rho_p = worker%get_field('Density(real)',    'value', 'face exterior')

        u_m   = worker%get_field('Velocity-1(real)', 'value', 'face interior')
        u_p   = worker%get_field('Velocity-1(real)', 'value', 'face exterior')

        v_m   = worker%get_field('Velocity-2(real)', 'value', 'face interior')
        v_p   = worker%get_field('Velocity-2(real)', 'value', 'face exterior')

        w_m   = worker%get_field('Velocity-3(real)', 'value', 'face interior')
        w_p   = worker%get_field('Velocity-3(real)', 'value', 'face exterior')

        p_m   = worker%get_field('Pressure(real)',   'value', 'face interior')
        p_p   = worker%get_field('Pressure(real)',   'value', 'face exterior')



        unorm_1 = worker%unit_normal_ale(1)
        unorm_2 = worker%unit_normal_ale(2)
        unorm_3 = worker%unit_normal_ale(3)




        !--------------------------------------
        !  Compute wave speeds
        !--------------------------------------


        !
        ! Compute normal velocities: dot-product vector projection along unit-normal direction
        !
        un = unorm_1*ubar + unorm_2*vbar + unorm_3*wbar

        !
        ! Compute wave speeds
        !
        wave_c = abs(un) + cbar


        wave = rho_m
        wave = wave_c


        !================================
        !       MASS FLUX
        !================================
        !upwind = -wave*(rho_p - rho_m)
        !integrand = HALF * ( upwind*normx*unormx  +  upwind*normy*unormy  +  upwind*normz*unormz )
        !call worker%integrate_boundary(irho,integrand)

        dissipation = HALF*wave*(rho_m - rho_p)

        call worker%integrate_boundary_upwind('Density(real)',dissipation)


        !================================
        !       X-MOMENTUM FLUX
        !================================
        !upwind = -wave*(u_p - u_m)
        !integrand = HALF * ( upwind*normx*unormx  +  upwind*normy*unormy  +  upwind*normz*unormz )
        !call worker%integrate_boundary(iu,integrand)

        dissipation = HALF*wave*(u_m - u_p)

        call worker%integrate_boundary_upwind('Velocity-1(real)',dissipation)

        !================================
        !       Y-MOMENTUM FLUX
        !================================
        !upwind = -wave*(v_p - v_m)
        !integrand = HALF * ( upwind*normx*unormx  +  upwind*normy*unormy  +  upwind*normz*unormz )
        !call worker%integrate_boundary(iv,integrand)

        dissipation = HALF*wave*(v_m - v_p)

        call worker%integrate_boundary_upwind('Velocity-2(real)',dissipation)

        !================================
        !       Z-MOMENTUM FLUX
        !================================
        !upwind = -wave*(w_p - w_m)
        !integrand = HALF * ( upwind*normx*unormx  +  upwind*normy*unormy  +  upwind*normz*unormz )
        !call worker%integrate_boundary(iw,integrand)

        dissipation = HALF*wave*(w_m - w_p)

        call worker%integrate_boundary_upwind('Velocity-3(real)',dissipation)


        !================================
        !          ENERGY FLUX
        !================================
        !upwind = -wave*(p_p - p_m)
        !integrand = HALF * ( upwind*normx*unormx  +  upwind*normy*unormy  +  upwind*normz*unormz )
        !call worker%integrate_boundary(ip,integrand)


        dissipation = HALF*wave*(w_m - w_p)

        call worker%integrate_boundary_upwind('Pressure(real)',dissipation)




    end subroutine compute
    !***************************************************************************************













end module PRIMLINEULER_LaxFriedrichs_flux_real
