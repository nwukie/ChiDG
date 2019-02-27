module HP_LaxFriedrichs
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF,THREE
    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/26/2019
    !!
    !-----------------------------------------------------------------------------------------------
    type, extends(operator_t), public :: HP_LaxFriedrichs_t

    contains

        procedure   :: init
        procedure   :: compute

    end type HP_LaxFriedrichs_t
    !***********************************************************************************************

    real(rk) :: p_param = 4._rk

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
            dissipation_u, dissipation_p, dissipation_q, dissipation_r, &
            sumsqr_m, sumsqr_p, mag_m, mag_p,   &
            d_norm_m, d_norm_p, lambda_m, lambda_p

        real(rk)    :: alpha, beta
            

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


        ! Compute normal for 1D system
        d_norm_m = p_m*worker%unit_normal(1) + q_m*worker%unit_normal(2) + r_m*worker%unit_normal(3)
        d_norm_p = p_p*worker%unit_normal(1) + q_p*worker%unit_normal(2) + r_p*worker%unit_normal(3)

        ! Eigenvalues for 1D p-Poisson system
        alpha = THREE + (p_param-TWO)/TWO
        beta  = TWO   + (p_param-TWO)/TWO

        lambda_m = sqrt(abs(alpha*d_norm_m**beta))
        lambda_p = sqrt(abs(alpha*d_norm_p**beta))




        !
        ! Compute dissipation
        !
        !dissipation_u = -HALF*max(mag_m,mag_p)*(u_p-u_m)
        !dissipation_p = -HALF*max(mag_m,mag_p)*(p_p-p_m)
        !dissipation_q = -HALF*max(mag_m,mag_p)*(q_p-q_m)
        !dissipation_r = -HALF*max(mag_m,mag_p)*(r_p-r_m)

        dissipation_u = -HALF*max(lambda_m,lambda_p)*(u_p-u_m)
        dissipation_p = -HALF*max(lambda_m,lambda_p)*(p_p-p_m)
        dissipation_q = -HALF*max(lambda_m,lambda_p)*(q_p-q_m)
        dissipation_r = -HALF*max(lambda_m,lambda_p)*(r_p-r_m)


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
