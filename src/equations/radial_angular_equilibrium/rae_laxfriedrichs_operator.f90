module rae_laxfriedrichs_operator
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ONE,TWO,HALF
    use mod_fluid,              only: omega
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none

    private


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/28/2016
    !!
    !!
    !----------------------------------------------------------------------------------
    type, extends(operator_t), public :: rae_laxfriedrichs_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rae_laxfriedrichs_operator_t
    !**********************************************************************************










contains

    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rae_laxfriedrichs_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name("RAE LaxFriedrichs Flux")

        ! Set operator type
        call self%set_operator_type("Boundary Advective Flux")

        ! Set operator equations
        call self%add_primary_field('Pressure-1')
        call self%add_primary_field('Pressure-2')

    end subroutine init
    !********************************************************************************









    !>  Compute Lax-Friedrichs upwind flux
    !!
    !!  Dissipation = -alpha(u_a_m - u_a_p)
    !!
    !!  Alpha is the maximum wave speed
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/16/2016
    !!
    !!------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(rae_laxfriedrichs_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::    &
            p1_m,   p1_p,   p2_m,   p2_p,           &
            wave,   upwind

        print*, 'laxfriedrichs - 1'

        !
        ! Interpolate solution to quadrature nodes
        !
        p1_m = worker%get_field('Pressure-1', 'value', 'face interior')
        p1_p = worker%get_field('Pressure-2', 'value', 'face exterior')

        p2_m = worker%get_field('Pressure-2', 'value', 'face interior')
        p2_p = worker%get_field('Pressure-2', 'value', 'face exterior')


        !
        ! Compute wave speeds
        !
        !wave = max(p1_m,p1_p,p2_m,p2_p)
        wave = 10._rk

        !===================================================
        !                   Momentum-1
        !===================================================
        !upwind = -HALF*wave*(p1_p - p1_m)
        upwind = HALF*wave*(p1_p - p1_m)

        call worker%integrate_boundary_upwind('Pressure-1',upwind)


        !===================================================
        !                   Momentum-2
        !===================================================
        !upwind = -HALF*wave*(p2_p - p2_m)
        upwind = HALF*wave*(p2_p - p2_m)

        call worker%integrate_boundary_upwind('Presure-2',upwind)


        print*, 'laxfriedrichs - 2'

    end subroutine compute
    !*******************************************************************************************




end module rae_laxfriedrichs_operator
