module GCL_LaxFriedrichs_operator
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
    type, extends(operator_t), public :: GCL_LaxFriedrichs_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type GCL_LaxFriedrichs_operator_t
    !***********************************************************************************************

contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine init(self)
        class(GCL_LaxFriedrichs_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('Geometric Conservation LaxFriedrichs Operator')

        ! Set operator type
        call self%set_operator_type('Boundary Advective Operator')

        ! Set operator equations
        call self%add_primary_field('g_bar')

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
        class(GCL_LaxFriedrichs_operator_t), intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        type(AD_D), dimension(:), allocatable   ::  &
            g_bar_m,  g_bar_p,                      &
            flux_1, flux_2, flux_3, integrand

        real(rk),   allocatable, dimension(:)   ::                      &
            grid_velocity_1, grid_velocity_2, grid_velocity_3,          &
            wave_speed, area, det_jacobian_grid_m, det_jacobian_grid_p, &
            norm_1, norm_2, norm_3, max_grid_vel

        real(rk),   allocatable, dimension(:,:,:)   ::  &
            inv_jacobian_grid_m, inv_jacobian_grid_p


        !
        ! Interpolate solution to quadrature nodes
        !
        g_bar_m = worker%get_primary_field_face('g_bar', 'value', 'face interior')
        g_bar_p = worker%get_primary_field_face('g_bar', 'value', 'face exterior')


        !
        ! Get model coefficients
        !
        grid_velocity_1     = worker%get_grid_velocity_face('u_grid','face interior')
        grid_velocity_2     = worker%get_grid_velocity_face('v_grid','face interior')
        grid_velocity_3     = worker%get_grid_velocity_face('w_grid','face interior')
        det_jacobian_grid_m = worker%get_det_jacobian_grid_face('value','face interior')
        det_jacobian_grid_p = worker%get_det_jacobian_grid_face('value','face exterior')
        inv_jacobian_grid_m = worker%get_inv_jacobian_grid_face('face interior')
        inv_jacobian_grid_p = worker%get_inv_jacobian_grid_face('face exterior')


        !
        ! Get normal vector
        !
        norm_1  = worker%normal(1)
        norm_2  = worker%normal(2)
        norm_3  = worker%normal(3)

               
        !
        ! Compute boundary upwind flux
        !
        max_grid_vel = max(abs(grid_velocity_1),abs(grid_velocity_2),abs(grid_velocity_3))
        wave_speed = max_grid_vel
        area = sqrt(norm_1**TWO+norm_2**TWO+norm_3**TWO)

        integrand = area*wave_speed*HALF*(g_bar_m-g_bar_p)


        !
        ! Integrate flux
        !
        integrand = -integrand
        call worker%integrate_boundary('g_bar',integrand)


    end subroutine compute
    !*************************************************************************************************








end module GCL_LaxFriedrichs_operator
