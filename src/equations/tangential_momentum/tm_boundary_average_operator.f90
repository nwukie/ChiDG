module tm_boundary_average_operator
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
    type, extends(operator_t), public :: tm_boundary_average_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type tm_boundary_average_operator_t
    !********************************************************************************










contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(tm_boundary_average_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("TM Boundary Average Flux")

        !
        ! Set operator type
        !
        call self%set_operator_type("Boundary Advective Flux")

        !
        ! Set operator equations
        !
        call self%add_primary_field("Pressure")

    end subroutine init
    !********************************************************************************



    !>  Boundary Flux routine for Euler
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!-------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(tm_boundary_average_operator_t), intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop

        ! Storage at quadrature nodes
        type(AD_D), allocatable,    dimension(:) :: &
            p_m,        p_p,                        &
            density_m,  density_p,                  &
            mom1_m,     mom1_p,                     &
            mom2_m,     mom2_p,                     &
            u_m,        u_p,                        &
            v_m,        v_p,                        &
            flux_1_m,   flux_2_m,   flux_3_m,       &
            flux_1_p,   flux_2_p,   flux_3_p

        real(rk),   allocatable, dimension(:)   :: r


        !
        ! Interpolate solution to quadrature nodes
        !
        p_m = worker%get_field('Pressure', 'value', 'face interior')
        p_p = worker%get_field('Pressure', 'value', 'face exterior')


        !
        ! Get model fields
        !
        density_m = worker%get_field('Density',    'value', 'face interior')
        mom1_m    = worker%get_field('Momentum-1', 'value', 'face interior')
        mom2_m    = worker%get_field('Momentum-2', 'value', 'face interior')

        density_p = worker%get_field('Density',    'value', 'face exterior')
        mom1_p    = worker%get_field('Momentum-1', 'value', 'face exterior')
        mom2_p    = worker%get_field('Momentum-2', 'value', 'face exterior')


        if (worker%coordinate_system() == 'Cylindrical') then
            r = worker%coordinate('1','face interior')
            mom2_m = mom2_m/r
            mom2_p = mom2_p/r
        end if




        u_m = mom1_m/density_m
        v_m = mom2_m/density_m

        u_p = mom1_p/density_p
        v_p = mom2_p/density_p


        !================================
        !           Momentum-1
        !================================
        flux_1_m = density_m
        flux_2_m = density_m*v_m*v_m  +  p_m
        flux_3_m = density_m
        flux_1_m = ZERO
        flux_3_m = ZERO

        flux_1_p = density_p
        flux_2_p = density_p*v_p*v_p  +  p_p
        flux_3_p = density_p
        flux_1_p = ZERO
        flux_3_p = ZERO

        call worker%integrate_boundary_average('Pressure','Advection',      &
                                                flux_1_m,flux_2_m,flux_3_m, &   
                                                flux_1_p,flux_2_p,flux_3_p)



    end subroutine compute
    !*********************************************************************************************************












end module tm_boundary_average_operator
