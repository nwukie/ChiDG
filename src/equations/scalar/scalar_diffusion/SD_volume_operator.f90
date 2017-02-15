module SD_volume_operator
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO,HALF

    use type_operator,          only: operator_t
    use type_chidg_worker,      only: chidg_worker_t
    use type_properties,        only: properties_t
    use DNAD_D

    use SD_properties,          only: SD_properties_t
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
    type, extends(operator_t), public :: SD_volume_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type SD_volume_operator_t
    !*************************************************************************

contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(SD_volume_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name("Scalar Diffusion Volume Operator")

        ! Set operator type
        call self%set_operator_type("Volume Diffusive Operator")

        ! Set operator equations
        call self%add_primary_field("u")

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
        class(SD_volume_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:)   ::  &
            flux_x, flux_y, flux_z, dudx, dudy, dudz, mu


        !
        ! Interpolate solution to quadrature nodes
        !
        dudx = worker%get_primary_field_element('u','grad1 + lift')
        dudy = worker%get_primary_field_element('u','grad2 + lift')
        dudz = worker%get_primary_field_element('u','grad3 + lift')


        !
        ! Compute scalar coefficient
        ! 
        mu = worker%get_model_field_element('Scalar Diffusion Coefficient', 'value')


        !
        ! Compute volume flux at quadrature nodes
        !
        flux_x = -mu*dudx
        flux_y = -mu*dudy
        flux_z = -mu*dudz


        !
        ! Integrate volume flux
        !
        call worker%integrate_volume('u',flux_x,flux_y,flux_z)



    end subroutine compute
    !****************************************************************************************************






end module SD_volume_operator
