module HP_boundary_average
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: ZERO,ONE,TWO,HALF
    use type_operator,              only: operator_t
    use type_chidg_worker,          only: chidg_worker_t
    use type_properties,            only: properties_t
    use DNAD_D
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: HP_boundary_average_t


    contains

        procedure   :: init
        procedure   :: compute

    end type HP_boundary_average_t
    !********************************************************************************


    real(rk) :: p_param = 4._rk

contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(HP_boundary_average_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('Hyperbolized Poisson Boundary Average Operator')

        ! Set operator type
        call self%set_operator_type('Boundary Advective Operator')

        ! Set operator equations
        call self%add_primary_field('u')
        call self%add_primary_field('p')
        call self%add_primary_field('q')
        call self%add_primary_field('r')

    end subroutine init
    !********************************************************************************










    !> Compute the average advective boundary flux for scalar linear advection
    !!
    !!   @author Nathan A. Wukie
    !!
    !!   @param[in]      mesh    Mesh data
    !!   @param[inout]   sdata   Solver data. Solution, RHS, Linearization etc.
    !!   @param[in]      ielem   Element index
    !!   @param[in]      iface   Face index
    !!   @param[in]      iblk    Block index indicating the linearization direction
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(HP_boundary_average_t),    intent(inout)   :: self
        type(chidg_worker_t),                               intent(inout)   :: worker
        class(properties_t),                                intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:)   ::  &
            u_m, u_p,                               &
            p_m, p_p,                               &
            q_m, q_p,                               &
            r_m, r_p,                               &
            flux_1_m, flux_2_m, flux_3_m,           &
            flux_1_p, flux_2_p, flux_3_p, sumsqr_m, sumsqr_p, mag_m, mag_p


        !
        ! Interpolate solution to quadrature nodes
        !
        u_m       = worker%get_field('u','value','face interior')
        u_p       = worker%get_field('u','value','face exterior')

        p_m = worker%get_field('p','value','face interior')
        p_p = worker%get_field('p','value','face exterior')
        q_m = worker%get_field('q','value','face interior')
        q_p = worker%get_field('q','value','face exterior')
        r_m = worker%get_field('r','value','face interior')
        r_p = worker%get_field('r','value','face exterior')


        ! Allocate derivatives
        flux_1_m = u_m*ZERO
        flux_2_m = u_m*ZERO
        flux_3_m = u_m*ZERO

        flux_1_p = u_p*ZERO
        flux_2_p = u_p*ZERO
        flux_3_p = u_p*ZERO


        !
        ! u-equation
        !

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

        flux_1_m = mag_m*p_m
        flux_2_m = mag_m*q_m
        flux_3_m = mag_m*r_m

        flux_1_p = mag_p*p_p
        flux_2_p = mag_p*q_p
        flux_3_p = mag_p*r_p


!        flux_1_m = p_m
!        flux_2_m = q_m
!        flux_3_m = r_m
!
!        flux_1_p = p_p
!        flux_2_p = q_p
!        flux_3_p = r_p

        call worker%integrate_boundary_average('u','Advection',                 &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)


        !
        ! p-equation
        !
        flux_1_m = u_m
        flux_2_m = ZERO
        flux_3_m = ZERO

        flux_1_p = u_p
        flux_2_p = ZERO
        flux_3_p = ZERO

        call worker%integrate_boundary_average('p','Advection',                 &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)

        !
        ! q-equation
        !
        flux_1_m = ZERO
        flux_2_m = u_m
        flux_3_m = ZERO

        flux_1_p = ZERO
        flux_2_p = u_p
        flux_3_p = ZERO

        call worker%integrate_boundary_average('q','Advection',                 &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)

        !
        ! r-equation
        !
        flux_1_m = ZERO
        flux_2_m = ZERO
        flux_3_m = u_m

        flux_1_p = ZERO
        flux_2_p = ZERO
        flux_3_p = u_p

        call worker%integrate_boundary_average('r','Advection',                 &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)

    end subroutine compute
    !**************************************************************************************************




end module HP_boundary_average
