module rac_laxfriedrichs_operator
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
    type, extends(operator_t), public :: rac_laxfriedrichs_operator_t

    contains

        procedure   :: init
        procedure   :: compute

    end type rac_laxfriedrichs_operator_t
    !**********************************************************************************










contains

    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(rac_laxfriedrichs_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name("RAC LaxFriedrichs Flux")

        ! Set operator type
        call self%set_operator_type("Boundary Advective Flux")

        ! Set operator equations
        call self%add_primary_field('Pressure')

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
        class(rac_laxfriedrichs_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop


        ! Storage at quadrature nodes
        type(AD_D), allocatable, dimension(:) ::    &
            p_m, p_p, wave, upwind


        !
        ! Interpolate solution to quadrature nodes
        !
        p_m = worker%get_field('Pressure', 'value', 'face interior')
        p_p = worker%get_field('Pressure', 'value', 'face exterior')

        !
        ! Compute wave speeds
        !
        wave = max(p_m,p_p)
        wave = 1.0_rk

        !===================================================
        !                   Momentum-1
        !===================================================
        upwind = -HALF*wave*(p_p - p_m)

        call worker%integrate_boundary_upwind('Pressure',upwind)



    end subroutine compute
    !*******************************************************************************************




end module rac_laxfriedrichs_operator
