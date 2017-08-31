module SA_LaxFriedrichs_operator
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
    type, extends(operator_t), public :: SA_LaxFriedrichs_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type SA_LaxFriedrichs_operator_t
    !***********************************************************************************************

contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine init(self)
        class(SA_LaxFriedrichs_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('Scalar Advection LaxFriedrichs Operator')

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
        class(SA_LaxFriedrichs_operator_t), intent(inout)   :: self
        type(chidg_worker_t),               intent(inout)   :: worker
        class(properties_t),                intent(inout)   :: prop


        type(AD_D), dimension(:), allocatable   ::  &
            u_m,  u_p,                              &
            c1_m, c2_m, c3_m,                       &
            c1_p, c2_p, c3_p,                       &
            c_n_m, c_n_p, dissipation, wave_speed
            
        real(rk),   allocatable, dimension(:,:) :: grid_vel
        real(rk),   allocatable, dimension(:)   :: unorm_1, unorm_2, unorm_3, grid_vel_n



        !
        ! Interpolate solution to quadrature nodes
        !
        u_m = worker%get_primary_field_face('u', 'value', 'face interior')
        u_p = worker%get_primary_field_face('u', 'value', 'face exterior')

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
        unorm_1  = worker%unit_normal_ale(1)
        unorm_2  = worker%unit_normal_ale(2)
        unorm_3  = worker%unit_normal_ale(3)

               
        !
        ! Compute boundary upwind dissipation
        !
        grid_vel   = worker%get_grid_velocity_face('face interior')
        c_n_m      = c1_m*unorm_1 + c2_m*unorm_2 + c3_m*unorm_3
        c_n_p      = c1_p*unorm_1 + c2_p*unorm_2 + c3_p*unorm_3
        grid_vel_n = grid_vel(:,1)*unorm_1  +  grid_vel(:,2)*unorm_2  +  grid_vel(:,3)*unorm_3

        wave_speed  = max(abs(c_n_m),abs(c_n_p)) + grid_vel_n
        dissipation = HALF*wave_speed*(u_m-u_p)


        !
        ! Integrate flux
        !
        call worker%integrate_boundary_upwind('u',dissipation)


    end subroutine compute
    !*************************************************************************************************








end module SA_LaxFriedrichs_operator
