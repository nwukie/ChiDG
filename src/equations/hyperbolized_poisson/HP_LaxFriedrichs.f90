module HP_LaxFriedrichs
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    type, extends(operator_t), public :: HP_LaxFriedrichs_t


    contains

        procedure   :: init
        procedure   :: compute

    end type HP_LaxFriedrichs_t
    !***********************************************************************************************

    real(rk) :: p_param = 2._rk

contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine init(self)
        class(HP_LaxFriedrichs_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('Hyperbolized Poisson LaxFriedrichs Operator')

        ! Set operator type
        call self%set_operator_type('Boundary Advective Operator')

        ! Set operator equations
        call self%add_primary_field('u')
        call self%add_primary_field('p')
        call self%add_primary_field('q')
        call self%add_primary_field('r')

    end subroutine init
    !***********************************************************************************************



    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(HP_LaxFriedrichs_t), intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        type(AD_D), dimension(:), allocatable   ::  &
            u_m, u_p,                               &
            p_m, p_p,                               &
            q_m, q_p,                               &
            r_m, r_p,                               &
            dissipation_u, dissipation_p, dissipation_q, dissipation_r, sumsqr_m, sumsqr_p, mag_m, mag_p
            

        !
        ! Interpolate solution to quadrature nodes
        !
        u_m = worker%get_field('u', 'value', 'face interior')
        u_p = worker%get_field('u', 'value', 'face exterior')
        
        p_m = worker%get_field('p', 'value', 'face interior')
        p_p = worker%get_field('p', 'value', 'face exterior')

        q_m = worker%get_field('q', 'value', 'face interior')
        q_p = worker%get_field('q', 'value', 'face exterior')

        r_m = worker%get_field('r', 'value', 'face interior')
        r_p = worker%get_field('r', 'value', 'face exterior')



        sumsqr_m = p_m*p_m + q_m*q_m + r_m*r_m
        sumsqr_p = p_p*p_p + q_p*q_p + r_p*r_p
        if (abs(p_param-2._rk) > 1.e-8_rk) then
            mag_m = sumsqr_m**((p_param-TWO)/TWO)
            mag_p = sumsqr_p**((p_param-TWO)/TWO)
        else
            mag_m = sumsqr_m
            mag_p = sumsqr_p
            mag_m = ONE
            mag_p = ONE
        end if




        !
        ! Compute dissipation
        !
        !wave_speed  = 1._rk
        !dissipation_u = -HALF*wave_speed*(u_p-u_m)
        !dissipation_p = -HALF*wave_speed*(p_p-p_m)
        dissipation_u = -HALF*max(mag_m,mag_p)*(u_p-u_m)
        dissipation_p = -HALF*max(mag_m,mag_p)*(p_p-p_m)
        dissipation_q = -HALF*max(mag_m,mag_p)*(q_p-q_m)
        dissipation_r = -HALF*max(mag_m,mag_p)*(r_p-r_m)



        !
        ! Integrate flux
        !
        call worker%integrate_boundary_upwind('u',dissipation_u)
        call worker%integrate_boundary_upwind('p',dissipation_p)
        call worker%integrate_boundary_upwind('q',dissipation_q)
        call worker%integrate_boundary_upwind('r',dissipation_r)


    end subroutine compute
    !*************************************************************************************************








end module HP_LaxFriedrichs
