module graddemo_P_volume_operator
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
    type, extends(operator_t), public :: graddemo_P_volume_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type graddemo_P_volume_operator_t
    !*************************************************************************

contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(graddemo_P_volume_operator_t),   intent(inout)      :: self

        ! Set operator name
        call self%set_name('Scalar Diffusion Volume Operator')

        ! Set operator type
        call self%set_operator_type('Volume Diffusive Operator')

        ! Set operator equations
        call self%add_primary_field('Pressure')

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
        class(graddemo_P_volume_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:)   ::  &
            grad1_pbc, grad2_pbc, grad3_pbc,        &
            grad1_p,   grad2_p,   grad3_p,          &
            flux_1, flux_2, flux_3
            


        !
        ! Interpolate solution to quadrature nodes
        !
        grad1_pbc = worker%get_field('Pressure', 'grad1', 'element')
        grad2_pbc = worker%get_field('Pressure', 'grad2', 'element')
        grad3_pbc = worker%get_field('Pressure', 'grad3', 'element')


!        grad1_p =
!        grad2_p =
!        grad3_p =


        !
        ! Compute volume flux at quadrature nodes
        !
        flux_1 = grad1_pbc - grad1_p
        flux_2 = grad2_pbc - grad2_p
        flux_3 = grad3_pbc - grad3_p


        !
        ! Integrate volume flux
        !
        call worker%integrate_volume_flux('Pressure','Diffusion',flux_1,flux_2,flux_3)



    end subroutine compute
    !****************************************************************************************************






end module graddemo_P_volume_operator
