module SD_ale_volume_operator
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D
    use ieee_arithmetic
    implicit none
    private

    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !!
    !-------------------------------------------------------------------------
    type, extends(operator_t), public :: SD_ale_volume_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type SD_ale_volume_operator_t
    !*************************************************************************

contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(SD_ale_volume_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name('Scalar Diffusion ALE Volume Operator')

        ! Set operator type
        call self%set_operator_type('Volume Diffusive Operator')

        ! Set operator equations
        call self%add_primary_field('u')

    end subroutine init
    !********************************************************************************





    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/14/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(SD_ale_volume_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:)   ::  &
            flux_1, flux_2, flux_3, mu

        type(AD_D), allocatable, dimension(:,:)   ::  &
            gradu, &
            flux_ref



        !
        ! Interpolate solution to quadrature nodes
        !
        gradu = worker%get_primary_field_grad_ale_element('u','gradient + lift')


        !
        ! Compute scalar coefficient
        ! 
        mu = worker%get_model_field_element('Scalar Diffusion Coefficient', 'value')


        !
        ! Compute volume flux at quadrature nodes
        !
        flux_1 = -mu*gradu(:,1)
        flux_2 = -mu*gradu(:,2)
        flux_3 = -mu*gradu(:,3)

        !
        ! Multiply by grid jacobian
        !
        flux_ref = worker%post_process_volume_diffusive_flux_ale(flux_1, flux_2, flux_3)

        !
        ! Integrate volume flux
        !
        call worker%integrate_volume('u',flux_ref(:,1),flux_ref(:,2),flux_ref(:,3))



    end subroutine compute
    !****************************************************************************************************






end module SD_ale_volume_operator
