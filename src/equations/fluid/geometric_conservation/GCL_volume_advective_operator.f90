module GCL_volume_advective_operator
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
    type, extends(operator_t), public :: GCL_volume_advective_operator_t


    contains
    
        procedure   :: init
        procedure   :: compute

    end type GCL_volume_advective_operator_t
    !*************************************************************************

contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(GCL_volume_advective_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('Geometric Conservation Volume Operator')

        ! Set operator type
        call self%set_operator_type('Volume Advective Operator')

        ! Set operator equations
        call self%add_primary_field('g_bar')

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
        class(GCL_volume_advective_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:)   ::  &
            g_bar, flux_1, flux_2, flux_3

        real(rk),   allocatable, dimension(:)     :: det_jacobian_grid
        real(rk),   allocatable, dimension(:,:)   :: grid_velocity
        real(rk),   allocatable, dimension(:,:,:) :: inv_jacobian_grid



        ! Just to set up AD variables
        g_bar = worker%get_field('g_bar','value','element')


        !
        ! Get model coefficients
        !
        grid_velocity     = worker%get_grid_velocity_element()
        det_jacobian_grid = worker%get_det_jacobian_grid_element('value')
        inv_jacobian_grid = worker%get_inv_jacobian_grid_element()



        !
        ! Compute volume flux at quadrature nodes
        !
        flux_1 = g_bar  !just to initialize AD allocation
        flux_2 = g_bar  !just to initialize AD allocation
        flux_3 = g_bar  !just to initialize AD allocation

        flux_1 = (inv_jacobian_grid(:,1,1)*grid_velocity(:,1) + inv_jacobian_grid(:,1,2)*grid_velocity(:,2) + inv_jacobian_grid(:,1,3)*grid_velocity(:,3))*det_jacobian_grid
        flux_2 = (inv_jacobian_grid(:,2,1)*grid_velocity(:,1) + inv_jacobian_grid(:,2,2)*grid_velocity(:,2) + inv_jacobian_grid(:,2,3)*grid_velocity(:,3))*det_jacobian_grid
        flux_3 = (inv_jacobian_grid(:,3,1)*grid_velocity(:,1) + inv_jacobian_grid(:,3,2)*grid_velocity(:,2) + inv_jacobian_grid(:,3,3)*grid_velocity(:,3))*det_jacobian_grid


        !
        ! Integrate volume flux
        !
        flux_1 = -flux_1
        flux_2 = -flux_2
        flux_3 = -flux_3
        call worker%integrate_volume_flux('g_bar','Advection',flux_1,flux_2,flux_3)

    end subroutine compute
    !****************************************************************************************************






end module GCL_volume_advective_operator
