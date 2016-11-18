module SA_boundary_average_advective_operator
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: ZERO,ONE,TWO,HALF, ME, NEIGHBOR

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

        integer(ik)                             :: iu

        type(AD_D), allocatable, dimension(:)   ::  &
            u_m, u_p,                               &
            dudx_m,dudy_m,dudz_m,                   &
            dudx_p,dudy_p,dudz_p,                   &
            cx_m, cy_m, cz_m,                       &
            cx_p, cy_p, cz_p,                       &
            flux_x, flux_y, flux_z, integrand

        real(rk),   allocatable, dimension(:)   ::  &
            normx, normy, normz


        !
        ! Get variable index
        !
        iu = prop%get_primary_field_index("u")

        
        !
        ! Interpolate solution to quadrature nodes
        !
        u_m    = worker%get_face_variable(iu, 'value' , ME)
        dudx_m = worker%get_face_variable(iu, 'ddx'   , ME)
        dudy_m = worker%get_face_variable(iu, 'ddy'   , ME)
        dudz_m = worker%get_face_variable(iu, 'ddz'   , ME)

        u_p    = worker%get_face_variable(iu, 'value' , NEIGHBOR)
        dudx_p = worker%get_face_variable(iu, 'ddx'   , NEIGHBOR)
        dudy_p = worker%get_face_variable(iu, 'ddy'   , NEIGHBOR)
        dudz_p = worker%get_face_variable(iu, 'ddz'   , NEIGHBOR)







        
        !
        ! Get model coefficients
        !
        !cx_m = prop%scalar%compute_cx(u_m)
        !cy_m = prop%scalar%compute_cy(u_m)
        !cz_m = prop%scalar%compute_cz(u_m)
        !cx_p = prop%scalar%compute_cx(u_p)
        !cy_p = prop%scalar%compute_cy(u_p)
        !cz_p = prop%scalar%compute_cz(u_p)
        cx_m = prop%scalar%compute_cx(u_m,dudx_m,dudy_m,dudz_m)
        cy_m = prop%scalar%compute_cy(u_m,dudx_m,dudy_m,dudz_m)
        cz_m = prop%scalar%compute_cz(u_m,dudx_m,dudy_m,dudz_m)
        cx_p = prop%scalar%compute_cx(u_p,dudx_p,dudy_p,dudz_p)
        cy_p = prop%scalar%compute_cy(u_p,dudx_p,dudy_p,dudz_p)
        cz_p = prop%scalar%compute_cz(u_p,dudx_p,dudy_p,dudz_p)


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
        call worker%integrate_boundary(iu,integrand)


    end subroutine compute
    !**************************************************************************************************




end module SA_boundary_average_advective_operator
