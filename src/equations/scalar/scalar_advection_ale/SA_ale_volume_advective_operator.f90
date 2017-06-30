module SA_ale_volume_advective_operator
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
    type, extends(operator_t), public :: SA_ale_volume_advective_operator_t


    contains
    
        procedure   :: init
        procedure   :: compute

    end type SA_ale_volume_advective_operator_t
    !*************************************************************************

contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(SA_ale_volume_advective_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('Scalar Advection ALE Volume Operator')

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
        class(SA_ale_volume_advective_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:)   ::  &
            u, flux_1, flux_2, flux_3, c1, c2, c3, flux_1_ref, flux_2_ref, flux_3_ref

        real(rk), allocatable, dimension(:) ::      &
           u_grid, v_grid, w_grid, det_jacobian_grid, testx


        real(rk), allocatable, dimension(:,:,:) ::      &
            jacobian_grid

        u_grid = worker%get_grid_velocity_element("u_grid")
        v_grid = worker%get_grid_velocity_element("v_grid")
        w_grid = worker%get_grid_velocity_element("w_grid")

!        print *, 'u_grid'
!        print *, u_grid
!        print *, 'v_grid'
!        print *, v_grid
!        print *, 'w_grid'
!        print *, w_grid
        jacobian_grid = worker%get_inv_jacobian_grid_element()
        det_jacobian_grid = worker%get_det_jacobian_grid_element()



        !
        ! Interpolate solution to quadrature nodes
        !
        u  = worker%get_primary_field_element('u','value')

        u = u/det_jacobian_grid

        !
        ! Get model coefficients
        !
        c1 = worker%get_model_field_element('Scalar Advection Velocity-1', 'value')
        c2 = worker%get_model_field_element('Scalar Advection Velocity-2', 'value')
        c3 = worker%get_model_field_element('Scalar Advection Velocity-3', 'value')


        !
        ! Compute volume flux at quadrature nodes
        !
        flux_1 = (c1 - u_grid ) * u
        flux_2 = (c2 - v_grid ) * u
        flux_3 = (c3 - w_grid ) * u

        flux_1_ref = det_jacobian_grid*(jacobian_grid(:,1,1)*flux_1 + jacobian_grid(:,1,2)*flux_2 + jacobian_grid(:,1,3)*flux_3)
        flux_2_ref = det_jacobian_grid*(jacobian_grid(:,2,1)*flux_1 + jacobian_grid(:,2,2)*flux_2 + jacobian_grid(:,2,3)*flux_3)
        flux_3_ref = det_jacobian_grid*(jacobian_grid(:,3,1)*flux_1 + jacobian_grid(:,3,2)*flux_2 + jacobian_grid(:,3,3)*flux_3)



        !
        ! Integrate volume flux
        !


        call worker%integrate_volume('u',flux_1_ref,flux_2_ref,flux_3_ref)

    end subroutine compute
    !****************************************************************************************************






end module SA_ale_volume_advective_operator
