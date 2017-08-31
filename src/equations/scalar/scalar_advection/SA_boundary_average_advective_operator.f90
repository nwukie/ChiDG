module SA_boundary_average_advective_operator
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: ZERO,ONE,TWO,HALF
    use type_operator,              only: operator_t
    use type_chidg_worker,          only: chidg_worker_t
    use type_properties,            only: properties_t
    use DNAD_D
    implicit none


    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: SA_boundary_average_advective_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type SA_boundary_average_advective_operator_t
    !********************************************************************************



contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(SA_boundary_average_advective_operator_t),   intent(inout)  :: self

        ! Set operator name
        call self%set_name('Scalar Advection Boundary Average Operator')

        ! Set operator type
        call self%set_operator_type('Boundary Advective Operator')

        ! Set operator equations
        call self%add_primary_field('u')

    end subroutine init
    !********************************************************************************










    !> Compute the average advective boundary flux for scalar linear advection
    !!
    !!   @author Nathan A. Wukie
    !!
    !!   @param[in]      mesh    Mesh data
    !!   @param[inout]   sdata   Solver data. Solution, RHS, Linearization etc.
    !!   @param[in]      ielem   Element index
    !!   @param[in]      iface   Face index
    !!   @param[in]      iblk    Block index indicating the linearization direction
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(SA_boundary_average_advective_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),                               intent(inout)   :: worker
        class(properties_t),                                intent(inout)   :: prop


        type(AD_D), allocatable, dimension(:)   ::  &
            u_m, u_p,                               &
            c1_m, c2_m, c3_m,                       &
            c1_p, c2_p, c3_p,                       &
            flux_1_m, flux_2_m, flux_3_m,           &
            flux_1_p, flux_2_p, flux_3_p

        
        !
        ! Interpolate solution to quadrature nodes
        !
        u_m = worker%get_primary_field_face('u','value' , 'face interior')
        u_p = worker%get_primary_field_face('u','value' , 'face exterior')


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
        ! Compute boundary average flux
        !
        flux_1_m = c1_m*u_m
        flux_2_m = c2_m*u_m
        flux_3_m = c3_m*u_m

        flux_1_p = c1_p*u_p
        flux_2_p = c2_p*u_p
        flux_3_p = c3_p*u_p

        !
        ! Integrate flux
        !
        call worker%integrate_boundary_average('u','Advection',                 &
                                                flux_1_m, flux_2_m, flux_3_m,   &
                                                flux_1_p, flux_2_p, flux_3_p)


    end subroutine compute
    !**************************************************************************************************




end module SA_boundary_average_advective_operator
