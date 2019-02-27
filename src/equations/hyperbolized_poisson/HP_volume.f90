module HP_volume
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
    !-------------------------------------------------------------------------
    type, extends(operator_t), public :: HP_volume_t


    contains
    
        procedure   :: init
        procedure   :: compute

    end type HP_volume_t
    !*************************************************************************

    real(rk) :: p_param = 4._rk

contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(HP_volume_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('Hyperbolized Poisson Volume Operator')

        ! Set operator type
        call self%set_operator_type('Volume Advective Operator')

        ! Set operator equations
        call self%add_primary_field('u')
        call self%add_primary_field('p')
        call self%add_primary_field('q')
        call self%add_primary_field('r')

    end subroutine init
    !********************************************************************************




    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(HP_volume_t),  intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:)   ::  &
            u, p, q, r, flux_1, flux_2, flux_3, sumsqr, mag


        !
        ! Interpolate solution to quadrature nodes
        !
        u  = worker%get_field('u','value','element')
        p  = worker%get_field('p','value','element')
        q  = worker%get_field('q','value','element')
        r  = worker%get_field('r','value','element')


        ! Allocate derivatives
        flux_1 = u*ZERO
        flux_2 = u*ZERO
        flux_3 = u*ZERO


        !
        ! u-equation
        !
        sumsqr = p*p + q*q + r*r
        if (abs(p_param-2._rk) > 1.e-8_rk) then
            mag = sumsqr**((p_param-TWO)/TWO)
        else
            mag = sumsqr
            mag = ONE
        end if

        flux_1 = mag*p
        flux_2 = mag*q
        flux_3 = mag*r

!        flux_1 = p
!        flux_2 = q
!        flux_3 = r

        call worker%integrate_volume_flux('u','Advection',flux_1,flux_2,flux_3)


        !
        ! p-equation
        !
        flux_1 = u
        flux_2 = ZERO
        flux_3 = ZERO

        call worker%integrate_volume_flux('p','Advection',flux_1,flux_2,flux_3)

        !
        ! q-equation
        !
        flux_1 = ZERO
        flux_2 = u
        flux_3 = ZERO

        call worker%integrate_volume_flux('q','Advection',flux_1,flux_2,flux_3)


        !
        ! r-equation
        !
        flux_1 = ZERO
        flux_2 = ZERO
        flux_3 = u

        call worker%integrate_volume_flux('r','Advection',flux_1,flux_2,flux_3)


    end subroutine compute
    !****************************************************************************************************






end module HP_volume
