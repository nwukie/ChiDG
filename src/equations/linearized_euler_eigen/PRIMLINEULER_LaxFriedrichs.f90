module PRIMLINEULER_LaxFriedrichs
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
    type, extends(operator_t), public :: PRIMLINEULER_LaxFriedrichs_t

    contains

        procedure   :: init
        procedure   :: compute

    end type PRIMLINEULER_LaxFriedrichs_t
    !**************************************************************************************




contains


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2018
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init(self)
        class(PRIMLINEULER_LaxFriedrichs_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('PRIMLINEULER LaxFriedrichs')

        ! Set operator type
        call self%set_operator_type('Boundary Advective Operator')

        ! Set operator equations
        call self%add_primary_field("Density"   )
        call self%add_primary_field("Velocity-1")
        call self%add_primary_field("Velocity-2")
        call self%add_primary_field("Velocity-3")
        call self%add_primary_field("Pressure"  )


    end subroutine init
    !*****************************************************************************************
    


    !>  Real component of numerical flux dissipation
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2018
    !!
    !--------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(PRIMLINEULER_LaxFriedrichs_t),    intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop



        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:)   ::  &
            rho_m,   rho_p,                         &
            u_m,     u_p,                           &
            v_m,     v_p,                           &
            w_m,     w_p,                           &
            p_m,     p_p,                           &
            dissipation,  wave, unorm_1, unorm_2, unorm_3, &
            un, wave_c


        !
        ! Zero primary matrix container to get just contribution from this integral
        !
        call worker%solverdata%lhs%clear()



        !
        ! Interpolate solution to quadrature nodes
        !
        rho_m = worker%get_field('Density',    'value', 'face interior')
        rho_p = worker%get_field('Density',    'value', 'face exterior')

        u_m   = worker%get_field('Velocity-1', 'value', 'face interior')
        u_p   = worker%get_field('Velocity-1', 'value', 'face exterior')

        v_m   = worker%get_field('Velocity-2', 'value', 'face interior')
        v_p   = worker%get_field('Velocity-2', 'value', 'face exterior')

        w_m   = worker%get_field('Velocity-3', 'value', 'face interior')
        w_p   = worker%get_field('Velocity-3', 'value', 'face exterior')

        p_m   = worker%get_field('Pressure',   'value', 'face interior')
        p_p   = worker%get_field('Pressure',   'value', 'face exterior')



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
        dissipation = HALF*wave*(rho_m - rho_p)

        call worker%integrate_boundary_upwind('Density',dissipation)


        !================================
        !       X-MOMENTUM FLUX
        !================================
        dissipation = HALF*wave*(u_m - u_p)

        call worker%integrate_boundary_upwind('Velocity-1',dissipation)

        !================================
        !       Y-MOMENTUM FLUX
        !================================
        dissipation = HALF*wave*(v_m - v_p)

        call worker%integrate_boundary_upwind('Velocity-2',dissipation)

        !================================
        !       Z-MOMENTUM FLUX
        !================================
        dissipation = HALF*wave*(w_m - w_p)

        call worker%integrate_boundary_upwind('Velocity-3',dissipation)


        !================================
        !          ENERGY FLUX
        !================================
        dissipation = HALF*wave*(w_m - w_p)

        call worker%integrate_boundary_upwind('Pressure',dissipation)


!        !
!        ! Store off matrix contributions
!        !
!        worker%solverdata%A2b = worker%solverdata%lhs



    end subroutine compute
    !***************************************************************************************













end module PRIMLINEULER_LaxFriedrichs
