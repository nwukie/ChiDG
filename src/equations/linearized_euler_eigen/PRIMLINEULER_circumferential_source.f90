module PRIMLINEULER_circumferential_source
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ONE,TWO,THREE,FOUR,FIVE,EIGHT,NINE,HALF,ZERO,PI

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
    type, extends(operator_t), public :: PRIMLINEULER_circumferential_source_t


    contains

        procedure   :: init
        procedure   :: compute

    end type PRIMLINEULER_circumferential_source_t
    !***********************************************************************************










contains

    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/22/2018
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(PRIMLINEULER_circumferential_source_t), intent(inout)   :: self

        ! Set operator name
        call self%set_name("PRIMLINEULER Circumferential Source")

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
        class(PRIMLINEULER_circumferential_source_t), intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop

        type(AD_D), allocatable, dimension(:)   ::  &
            rho_r, u_r, v_r, w_r, p_r,              &
            source

        real(rk),   allocatable, dimension(:)   :: r

        print*, 'Computing C-source'

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

        !
        ! Get radial coordinates
        !
        r = worker%coordinate('1','volume')


        !===========================
        !        MASS FLUX
        !===========================
        source = (rho_2_rho * rho_r  + &
                  rho_2_u   * u_r    + &
                  rho_2_v   * v_r    + &
                  rho_2_w   * w_r    + &
                  rho_2_p   * p_r)*(mod_m/r)


        call worker%integrate_volume_source('Density',source)



        !===========================
        !     X-MOMENTUM FLUX
        !===========================
        source = (u_2_rho * rho_r  + &
                  u_2_u   * u_r    + &
                  u_2_v   * v_r    + &
                  u_2_w   * w_r    + &
                  u_2_p   * p_r)*(mod_m/r)

        call worker%integrate_volume_source('Velocity-1',source)


        !============================
        !     Y-MOMENTUM FLUX
        !============================
        source = (v_2_rho * rho_r  + &
                  v_2_u   * u_r    + &
                  v_2_v   * v_r    + &
                  v_2_w   * w_r    + &
                  v_2_p   * p_r)*(mod_m/r)

        call worker%integrate_volume_source('Velocity-2',source)

        !============================
        !     Z-MOMENTUM FLUX
        !============================
        source = (w_2_rho * rho_r  + &
                  w_2_u   * u_r    + &
                  w_2_v   * v_r    + &
                  w_2_w   * w_r    + &
                  w_2_p   * p_r)*(mod_m/r)

        call worker%integrate_volume_source('Velocity-3',source)

        !============================
        !       ENERGY FLUX
        !============================
        source = (p_2_rho * rho_r  + &
                  p_2_u   * u_r    + &
                  p_2_v   * v_r    + &
                  p_2_w   * w_r    + &
                  p_2_p   * p_r)*(mod_m/r)

        call worker%integrate_volume_source('Pressure',source)



!        !
!        ! Store off matrix contributions
!        !
!        worker%solverdata%B = worker%solverdata%lhs

    end subroutine compute
    !******************************************************************************************************






end module PRIMLINEULER_circumferential_source
