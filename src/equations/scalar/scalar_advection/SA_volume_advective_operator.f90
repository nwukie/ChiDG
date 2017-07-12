module SA_volume_advective_operator
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
    type, extends(operator_t), public :: SA_volume_advective_operator_t


    contains
    
        procedure   :: init
        procedure   :: compute

    end type SA_volume_advective_operator_t
    !*************************************************************************

contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(SA_volume_advective_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('Scalar Advection Volume Operator')

        ! Set operator type
        call self%set_operator_type('Volume Advective Operator')

        ! Set operator equations
        call self%add_primary_field('u')

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
        class(SA_volume_advective_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:)   ::  &
            u, flux_1, flux_2, flux_3, c1, c2, c3


        !
        ! Interpolate solution to quadrature nodes
        !
        !u  = worker%get_primary_field_element('u','value')
        u  = worker%get_field('u','value','element')


        !
        ! Get model coefficients
        !
        !c1 = worker%get_model_field_element('Scalar Advection Velocity-1', 'value')
        !c2 = worker%get_model_field_element('Scalar Advection Velocity-2', 'value')
        !c3 = worker%get_model_field_element('Scalar Advection Velocity-3', 'value')
        c1 = worker%get_field('Scalar Advection Velocity-1', 'value', 'element')
        c2 = worker%get_field('Scalar Advection Velocity-2', 'value', 'element')
        c3 = worker%get_field('Scalar Advection Velocity-3', 'value', 'element')


        !
        ! Compute volume flux at quadrature nodes
        !
        flux_1 = c1 * u
        flux_2 = c2 * u
        flux_3 = c3 * u


        !
        ! Integrate volume flux
        !
        call worker%integrate_volume('u',flux_1,flux_2,flux_3)



    end subroutine compute
    !****************************************************************************************************






end module SA_volume_advective_operator
