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
        call self%set_name("Scalar Advection Boundary Average Operator")

        ! Set operator type
        call self%set_operator_type("Boundary Advective Operator")

        ! Set operator equations
        call self%add_primary_field("u")

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
            cx_m, cy_m, cz_m,                       &
            cx_p, cy_p, cz_p,                       &
            flux_x, flux_y, flux_z, integrand

        real(rk),   allocatable, dimension(:)   ::  &
            normx, normy, normz


        
        !
        ! Interpolate solution to quadrature nodes
        !
        u_m = worker%get_primary_field_face('u','value' , 'face interior')
        u_p = worker%get_primary_field_face('u','value' , 'face exterior')

        
        !
        ! Get model coefficients
        !
        cx_m = worker%get_model_field_face('Scalar X-Advection Velocity', 'value', 'face interior')
        cy_m = worker%get_model_field_face('Scalar Y-Advection Velocity', 'value', 'face interior')
        cz_m = worker%get_model_field_face('Scalar Z-Advection Velocity', 'value', 'face interior')
        cx_p = worker%get_model_field_face('Scalar X-Advection Velocity', 'value', 'face exterior')
        cy_p = worker%get_model_field_face('Scalar Y-Advection Velocity', 'value', 'face exterior')
        cz_p = worker%get_model_field_face('Scalar Z-Advection Velocity', 'value', 'face exterior')


        !
        ! Get normal vector
        !
        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)


        !
        ! Compute boundary average flux
        !
        flux_x = HALF*(cx_m*u_m + cx_p*u_p)
        flux_y = HALF*(cy_m*u_m + cy_p*u_p)
        flux_z = HALF*(cz_m*u_m + cz_p*u_p)


        !
        ! Dot with normal vector
        ! 
        integrand = flux_x*normx + flux_y*normy + flux_z*normz


        !
        ! Integrate flux
        !
        call worker%integrate_boundary('u',integrand)


    end subroutine compute
    !**************************************************************************************************




end module SA_boundary_average_advective_operator
