module SA_ale_LaxFriedrichs_operator
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
    type, extends(operator_t), public :: SA_ale_LaxFriedrichs_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type SA_ale_LaxFriedrichs_operator_t
    !***********************************************************************************************

contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine init(self)
        class(SA_ale_LaxFriedrichs_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('Scalar Advection ALE LaxFriedrichs Operator')

        ! Set operator type
        call self%set_operator_type('Boundary Advective Operator')

        ! Set operator equations
        call self%add_primary_field('u')

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
        class(SA_ale_LaxFriedrichs_operator_t), intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        type(AD_D), dimension(:), allocatable   ::  &
            u_m,  u_p,                              &
            c1_m, c2_m, c3_m,                       &
            c1_p, c2_p, c3_p,                       &
            flux_1, flux_2, flux_3, integrand

        type(AD_D), allocatable, dimension(:,:)   :: flux_ref

        real(rk), allocatable                       :: u_grid(:), v_grid(:), w_grid(:),wave_speed(:), area(:)
        
        real(rk),   allocatable, dimension(:)   ::  &
            norm_1, norm_2, norm_3


        u_grid = worker%get_grid_velocity_face('u_grid','face interior')
        v_grid = worker%get_grid_velocity_face('v_grid','face interior')
        w_grid = worker%get_grid_velocity_face('w_grid','face interior')

        !
        ! Interpolate solution to quadrature nodes
        !
        u_m = worker%get_primary_field_value_ale_face('u', 'face interior')
        u_p = worker%get_primary_field_value_ale_face('u', 'face exterior')

        !
        ! Get model coefficients
        !
        c1_m = worker%get_model_field_face('Scalar Advection Velocity-1', 'value', 'face interior')
        c2_m = worker%get_model_field_face('Scalar Advection Velocity-2', 'value', 'face interior')
        c3_m = worker%get_model_field_face('Scalar Advection Velocity-3', 'value', 'face interior')
        c1_p = worker%get_model_field_face('Scalar Advection Velocity-1', 'value', 'face exterior')
        c2_p = worker%get_model_field_face('Scalar Advection Velocity-2', 'value', 'face exterior')
        c3_p = worker%get_model_field_face('Scalar Advection Velocity-3', 'value', 'face exterior')
        

        !
        ! Get normal vector
        !
        norm_1  = worker%normal(1)
        norm_2  = worker%normal(2)
        norm_3  = worker%normal(3)

               
        !
        ! Compute boundary upwind flux
        !
        !flux_1 = max(abs(c1_m)+abs(u_grid),abs(c1_p)+abs(u_grid))*HALF*(u_m - u_p)
        !flux_2 = max(abs(c2_m)+abs(v_grid),abs(c2_p)+abs(v_grid))*HALF*(u_m - u_p)
        !flux_3 = max(abs(c3_m)+abs(w_grid),abs(c3_p)+abs(w_grid))*HALF*(u_m - u_p)

        !flux_ref = worker%post_process_boundary_diffusive_flux_ale(flux_1, flux_2, flux_3)
        
        !
        ! Dot with normal vector
        ! 

        !integrand = flux_ref(:,1)*norm_1 + flux_ref(:,2)*norm_2 + flux_ref(:,3)*norm_3
        !max_grid_vel = max(abs(u_grid),abs(v_grid),abs(w_grid))
        wave_speed = ONE + abs(u_grid)
        area = sqrt(norm_1**TWO+norm_2**TWO+norm_3**TWO)
        integrand = area*wave_speed*HALF*(u_m-u_p)


        !
        ! Integrate flux
        !
        call worker%integrate_boundary('u',integrand)


    end subroutine compute
    !*************************************************************************************************








end module SA_ale_LaxFriedrichs_operator
