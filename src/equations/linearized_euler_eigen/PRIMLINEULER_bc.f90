module PRIMLINEULER_bc
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ONE,TWO,THREE,FOUR,FIVE,EIGHT,NINE,HALF,ZERO,PI

    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D

    use mod_primitive_linearized_euler
    implicit none
    private


    !>  BC flux for Linearized Euler equations.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2018
    !!
    !!
    !----------------------------------------------------------------------------------
    type, extends(operator_t), public :: PRIMLINEULER_bc_t

    contains

        procedure   :: init
        procedure   :: compute

    end type PRIMLINEULER_bc_t
    !***********************************************************************************










contains

    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(PRIMLINEULER_bc_t),  intent(inout)   :: self

        ! Set operator name
        call self%set_name("PRIMLINEULER BC")

        ! Set operator type
        call self%set_operator_type("BC Advective Flux")

        ! Set operator equations
        call self%add_primary_field("Density"   )
        call self%add_primary_field("Velocity-1")
        call self%add_primary_field("Velocity-2")
        call self%add_primary_field("Velocity-3")
        call self%add_primary_field("Pressure"  )

    end subroutine init
    !********************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/23/2018
    !!
    !----------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(PRIMLINEULER_bc_t),   intent(inout)   :: self
        type(chidg_worker_t),       intent(inout)   :: worker
        class(properties_t),        intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:)   ::  &
            rho_r, u_r, v_r, w_r, p_r,              &
            flux_1, flux_2, flux_3

        print*, 'Computing BC'

        !
        ! Zero primary matrix container to get just contribution from this integral
        !
        call worker%solverdata%lhs%clear()




        !
        ! Interpolate solution to quadrature nodes: REAL
        !
        rho_r = worker%get_field('Density'   , 'value', 'face interior')
        u_r   = worker%get_field('Velocity-1', 'value', 'face interior')
        v_r   = worker%get_field('Velocity-2', 'value', 'face interior')
        w_r   = worker%get_field('Velocity-3', 'value', 'face interior')
        p_r   = worker%get_field('Pressure'  , 'value', 'face interior')



        !===========================
        !        MASS FLUX
        !===========================
        flux_1 = rho_1_rho * rho_r  + &
!                 rho_1_u   * u_r    + &
!                 rho_1_v   * v_r    + &
!                 rho_1_w   * w_r    + &
                 rho_1_p   * p_r

        flux_2 = rho_2_rho * rho_r  + &
!                 rho_2_u   * u_r    + &
!                 rho_2_v   * v_r    + &
!                 rho_2_w   * w_r    + &
                 rho_2_p   * p_r

        flux_3 = rho_3_rho * rho_r  + &
!                 rho_3_u   * u_r    + &
!                 rho_3_v   * v_r    + &
!                 rho_3_w   * w_r    + &
                 rho_3_p   * p_r

        flux_2 = ZERO
        flux_3 = ZERO
        call worker%integrate_boundary_condition('Density','Advection',flux_1,flux_2,flux_3)



        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        flux_1 = u_1_rho * rho_r  + &
!                 u_1_u   * u_r    + &
!                 u_1_v   * v_r    + &
!                 u_1_w   * w_r    + &
                 u_1_p   * p_r

        flux_2 = u_2_rho * rho_r  + &
!                 u_2_u   * u_r    + &
!                 u_2_v   * v_r    + &
!                 u_2_w   * w_r    + &
                 u_2_p   * p_r

        flux_3 = u_3_rho * rho_r  + &
!                 u_3_u   * u_r    + &
!                 u_3_v   * v_r    + &
!                 u_3_w   * w_r    + &
                 u_3_p   * p_r

        flux_2 = ZERO
        flux_3 = ZERO
        call worker%integrate_boundary_condition('Velocity-1','Advection',flux_1,flux_2,flux_3)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        flux_1 = v_1_rho * rho_r  + &
!                 v_1_u   * u_r    + &
!                 v_1_v   * v_r    + &
!                 v_1_w   * w_r    + &
                 v_1_p   * p_r

        flux_2 = v_2_rho * rho_r  + &
!                 v_2_u   * u_r    + &
!                 v_2_v   * v_r    + &
!                 v_2_w   * w_r    + &
                 v_2_p   * p_r

        flux_3 = v_3_rho * rho_r  + &
!                 v_3_u   * u_r    + &
!                 v_3_v   * v_r    + &
!                 v_3_w   * w_r    + &
                 v_3_p   * p_r

        flux_2 = ZERO
        flux_3 = ZERO
        call worker%integrate_boundary_condition('Velocity-2','Advection',flux_1,flux_2,flux_3)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        flux_1 = w_1_rho * rho_r  + &
!                 w_1_u   * u_r    + &
!                 w_1_v   * v_r    + &
!                 w_1_w   * w_r    + &
                 w_1_p   * p_r

        flux_2 = w_2_rho * rho_r  + &
!                 w_2_u   * u_r    + &
!                 w_2_v   * v_r    + &
!                 w_2_w   * w_r    + &
                 w_2_p   * p_r

        flux_3 = w_3_rho * rho_r  + &
!                 w_3_u   * u_r    + &
!                 w_3_v   * v_r    + &
!                 w_3_w   * w_r    + &
                 w_3_p   * p_r

        flux_2 = ZERO
        flux_3 = ZERO
        call worker%integrate_boundary_condition('Velocity-3','Advection',flux_1,flux_2,flux_3)

        !============================
        !       ENERGY FLUX
        !============================
        flux_1 = p_1_rho * rho_r  + &
!                 p_1_u   * u_r    + &
!                 p_1_v   * v_r    + &
!                 p_1_w   * w_r    + &
                 p_1_p   * p_r

        flux_2 = p_2_rho * rho_r  + &
!                 p_2_u   * u_r    + &
!                 p_2_v   * v_r    + &
!                 p_2_w   * w_r    + &
                 p_2_p   * p_r

        flux_3 = p_3_rho * rho_r  + &
!                 p_3_u   * u_r    + &
!                 p_3_v   * v_r    + &
!                 p_3_w   * w_r    + &
                 p_3_p   * p_r

        flux_2 = ZERO
        flux_3 = ZERO
        call worker%integrate_boundary_condition('Pressure','Advection',flux_1,flux_2,flux_3)


!        !
!        ! Store off matrix contributions
!        !
!        if (worker%iface == 1) then
!            worker%solverdata%A3a = worker%solverdata%lhs
!        else if (worker%iface == 2) then
!            worker%solverdata%A3b = worker%solverdata%lhs
!        end if



    end subroutine compute
    !******************************************************************************************************






end module PRIMLINEULER_bc
