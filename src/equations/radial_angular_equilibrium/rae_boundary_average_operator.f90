module rae_boundary_average_operator
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE, TWO, HALF, ZERO
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none

    private



    !> Implementation of the Euler boundary average flux
    !!
    !!  - At a boundary interface, the solution states Q- and Q+ exists on opposite 
    !!    sides of the boundary. The average flux is computed as Favg = 1/2(F(Q-) + F(Q+))
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: rae_boundary_average_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rae_boundary_average_operator_t
    !********************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rae_boundary_average_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("RAE Boundary Average Flux")

        !
        ! Set operator type
        !
        call self%set_operator_type("Boundary Advective Flux")

        !
        ! Set operator equations
        !
        call self%add_primary_field("Pressure-1")
        call self%add_primary_field("Pressure-2")

    end subroutine init
    !********************************************************************************



    !>  Boundary Flux routine for Euler
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!-------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(rae_boundary_average_operator_t), intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable,    dimension(:) :: &
            p1_m,       p1_p,                       &
            p2_m,       p2_p,                       &
            density_m,  density_p,                  &
            u_m,        u_p,                        &
            v_m,        v_p,                        &
            flux_1_m,   flux_2_m,   flux_3_m,       &
            flux_1_p,   flux_2_p,   flux_3_p

        print*, 'bounday average - 1'

        !
        ! Interpolate solution to quadrature nodes
        !
        p1_m = worker%get_field('Pressure-1', 'value', 'face interior')
        p1_p = worker%get_field('Pressure-1', 'value', 'face exterior')

        p2_m = worker%get_field('Pressure-2', 'value', 'face interior')
        p2_p = worker%get_field('Pressure-2', 'value', 'face exterior')


        !
        ! Get model fields
        !
        density_m = worker%get_field('Density',    'value', 'face interior')
        u_m       = worker%get_field('Velocity-1', 'value', 'face interior')
        v_m       = worker%get_field('Velocity-2', 'value', 'face interior')

        density_p = worker%get_field('Density',    'value', 'face exterior')
        u_p       = worker%get_field('Velocity-1', 'value', 'face exterior')
        v_p       = worker%get_field('Velocity-2', 'value', 'face exterior')



        !================================
        !           Momentum-1
        !================================
        flux_1_m = (density_m*u_m * u_m)  +  p1_m*p2_m
        flux_2_m = (density_m*u_m * v_m)
        flux_3_m = (density_m*u_m)
        flux_3_m = ZERO

        flux_1_p = (density_p*u_p * u_p)  +  p1_p*p2_p
        flux_2_p = (density_p*u_p * v_p)
        flux_3_p = (density_p*u_p)
        flux_3_p = ZERO

        call worker%integrate_boundary_average('Pressure-1','Advection',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)


        !================================
        !           Momentum-2 
        !================================
        flux_1_m = (density_m*v_m * u_m)
        flux_2_m = (density_m*v_m * v_m)  +  p1_m*p2_m
        flux_3_m = (density_m*v_m)
        flux_3_m = ZERO

        flux_1_p = (density_p*v_p * u_p)
        flux_2_p = (density_p*v_p * v_p)  +  p1_p*p2_p
        flux_3_p = (density_p*v_p)
        flux_3_p = ZERO

        call worker%integrate_boundary_average('Pressure-2','Advection',    &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)

        print*, 'bounday average - 2'

    end subroutine compute
    !*********************************************************************************************************












end module rae_boundary_average_operator
