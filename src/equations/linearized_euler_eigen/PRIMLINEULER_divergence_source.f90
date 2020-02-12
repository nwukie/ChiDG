module PRIMLINEULER_divergence_source
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ONE,ZERO

    use type_operator,      only: operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    use DNAD_D

    use mod_primitive_linearized_euler
    implicit none
    private


    !>  Volume advective flux for Linearized Euler equations - real.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2018
    !!
    !!
    !----------------------------------------------------------------------------------
    type, extends(operator_t), public :: PRIMLINEULER_divergence_source_t


    contains

        procedure   :: init
        procedure   :: compute

    end type PRIMLINEULER_divergence_source_t
    !***********************************************************************************










contains

    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(PRIMLINEULER_divergence_source_t), intent(inout)   :: self

        ! Set operator name
        call self%set_name("PRIMLINEULER Divergence Source")

        ! Set operator type
        call self%set_operator_type("Volume Advective Flux")

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
        class(PRIMLINEULER_divergence_source_t),    intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:)   ::  &
            rho_r, u_r, v_r, w_r, p_r, source, r


        print*, 'Computing Divergence Source'

        !
        ! Zero primary matrix container to get just contribution from this integral
        !
        call worker%solverdata%lhs%clear()



        !
        ! Interpolate solution to quadrature nodes: REAL
        !
        rho_r = worker%get_field('Density'   , 'value', 'element')
        u_r   = worker%get_field('Velocity-1', 'value', 'element')
        v_r   = worker%get_field('Velocity-2', 'value', 'element')
        w_r   = worker%get_field('Velocity-3', 'value', 'element')
        p_r   = worker%get_field('Pressure'  , 'value', 'element')


        r = worker%coordinate('1', 'volume')


        !===========================
        !        MASS FLUX
        !===========================
        source = (ONE/r)*(rho_1_rho * rho_r  + &
                          rho_1_u   * u_r    + &
                          rho_1_v   * v_r    + &
                          rho_1_w   * w_r    + &
                          rho_1_p   * p_r)

        call worker%integrate_volume_source('Density',source)



        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        source = (ONE/r)*(u_1_rho * rho_r  + &
                          u_1_u   * u_r    + &
                          u_1_v   * v_r    + &
                          u_1_w   * w_r    + &
                          u_1_p   * p_r)
        call worker%integrate_volume_source('Velocity-1',source)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        source = (ONE/r)*(v_1_rho * rho_r  + &
                          v_1_u   * u_r    + &
                          v_1_v   * v_r    + &
                          v_1_w   * w_r    + &
                          v_1_p   * p_r)

        call worker%integrate_volume_source('Velocity-2',source)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        source = (ONE/r)*(w_1_rho * rho_r  + &
                          w_1_u   * u_r    + &
                          w_1_v   * v_r    + &
                          w_1_w   * w_r    + &
                          w_1_p   * p_r)

        call worker%integrate_volume_source('Velocity-3',source)

        !============================
        !       ENERGY FLUX
        !============================
        source = (ONE/r)*(p_1_rho * rho_r  + &
                          p_1_u   * u_r    + &
                          p_1_v   * v_r    + &
                          p_1_w   * w_r    + &
                          p_1_p   * p_r)

        call worker%integrate_volume_source('Pressure',source)




!        !
!        ! Store off matrix contributions
!        !
!        worker%solverdata%F = worker%solverdata%lhs

    end subroutine compute
    !******************************************************************************************************






end module PRIMLINEULER_divergence_source
