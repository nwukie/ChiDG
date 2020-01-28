module DLA_LaxFriedrichs_flux
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ZERO,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG, &
                                      ME, NEIGHBOR

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none


    !> This equation set exists really just to test equationsets with 
    !! more than one equation. The idea is just to compute the linear
    !! advection solution twice at the same time. The equations are 
    !! independent of each other. So, we can verify, for example,
    !! the volume flux jacobians for each equation. They should be the
    !! same as for the single LinearAdvection equation set
    !!
    !!  @author Nathan A. Wukie
    !!
    !-----------------------------------------------------------------------------------------------
    type, extends(operator_t), public :: DLA_LaxFriedrichs_flux_t

    contains

        procedure   :: init
        procedure   :: compute

    end type DLA_LaxFriedrichs_flux_t
    !************************************************************************************************

contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(DLA_LaxFriedrichs_flux_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name("DLA LaxFriedrichs Flux")

        ! Set operator type
        call self%set_operator_type("Boundary Advective Flux")

        ! Set operator equations
        call self%add_primary_field("u_a")
        call self%add_primary_field("u_b")

    end subroutine init
    !********************************************************************************








    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(DLA_LaxFriedrichs_flux_t),    intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:)   ::  & 
            ua_p, ua_m, ub_p, ub_m, dissipation, wave_speed, c1, c2, c3, c_n_m, c_n_p,  &
            unorm_1, unorm_2, unorm_3, grid_vel_n

        real(rk),   allocatable, dimension(:,:) :: grid_vel


        ! Interpolate solution to quadrature nodes
        ua_m = worker%get_field('u_a','value', 'face interior')
        ua_p = worker%get_field('u_a','value', 'face exterior')

        ub_m = worker%get_field('u_b','value', 'face interior')
        ub_p = worker%get_field('u_b','value', 'face exterior')

        ! Get equation set property data
        c1 = ua_m
        c2 = ua_m
        c3 = ua_m
        c1 = ONE
        c2 = ZERO
        c3 = ZERO

        ! Get normal vector
        unorm_1  = worker%unit_normal_ale(1)
        unorm_2  = worker%unit_normal_ale(2)
        unorm_3  = worker%unit_normal_ale(3)

               
        ! Compute boundary upwind dissipation
        grid_vel   = worker%get_grid_velocity_face('face interior')
        c_n_m      = c1*unorm_1 + c2*unorm_2 + c3*unorm_3
        c_n_p      = c1*unorm_1 + c2*unorm_2 + c3*unorm_3
        grid_vel_n = grid_vel(:,1)*unorm_1  +  grid_vel(:,2)*unorm_2  +  grid_vel(:,3)*unorm_3
        wave_speed  = max(abs(c_n_m),abs(c_n_p)) + grid_vel_n


        ! Integrate flux
        dissipation = HALF*wave_speed*(ua_m-ua_p)
        call worker%integrate_boundary_upwind('u_a',dissipation)

        dissipation = HALF*wave_speed*(ub_m-ub_p)
        call worker%integrate_boundary_upwind('u_b',dissipation)


    end subroutine compute
    !****************************************************************************************








end module DLA_LaxFriedrichs_flux
