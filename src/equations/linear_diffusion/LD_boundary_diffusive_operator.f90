module LD_boundary_diffusive_operator
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: ZERO,ONE,TWO,HALF, ME, NEIGHBOR

    use type_operator,              only: operator_t
    use type_chidg_worker,          only: chidg_worker_t
    use type_properties,            only: properties_t
    use DNAD_D

    use LD_properties,              only: LD_properties_t
    implicit none

    private



    !>
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------
    type, extends(operator_t), public :: LD_boundary_diffusive_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type LD_boundary_diffusive_operator_t
    !********************************************************************************



contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(LD_boundary_diffusive_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Linear Diffusion Boundary Average Flux")

        !
        ! Set operator type
        !
        call self%set_operator_type("Boundary Diffusive Flux")

        !
        ! Set operator equations
        !
        call self%set_equation("u")

    end subroutine init
    !********************************************************************************







    !>  Compute the diffusive boundary flux for scalar linear diffusion.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      mesh    Mesh data
    !!  @param[inout]   sdata   Solver data. Solution, RHS, Linearization etc.
    !!  @param[in]      ielem   Element index
    !!  @param[in]      iface   Face index
    !!  @param[in]      iblk    Block index indicating the linearization direction
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute(self,worker,prop)
        class(LD_boundary_diffusive_operator_t),    intent(inout)   :: self
        type(chidg_worker_t),                       intent(inout)   :: worker
        class(properties_t),                        intent(inout)   :: prop


        integer(ik)                             :: iu

        type(AD_D), allocatable, dimension(:)   :: &
            dudx_m, dudy_m, dudz_m, dudx_p, dudy_p, dudz_p, &
            flux_x, flux_y, flux_z, flux_m, flux_p, flux, integrand

        real(rk),   allocatable, dimension(:)   :: &
            normx, normy, normz


        !
        ! Get variable index
        !
        iu = prop%get_equation_index("u")

        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)


        !
        ! Interpolate solution to quadrature nodes
        !
        dudx_m = worker%get_face_variable(iu, 'ddx + lift', ME)
        dudy_m = worker%get_face_variable(iu, 'ddy + lift', ME)
        dudz_m = worker%get_face_variable(iu, 'ddz + lift', ME)


        dudx_p = worker%get_face_variable(iu, 'ddx + lift', NEIGHBOR)
        dudy_p = worker%get_face_variable(iu, 'ddy + lift', NEIGHBOR)
        dudz_p = worker%get_face_variable(iu, 'ddz + lift', NEIGHBOR)


        flux_m = -dudx_m
        flux_p = -dudx_p
        flux_x = HALF*(flux_m + flux_p)

        flux_m = -dudy_m
        flux_p = -dudy_p
        flux_y = HALF*(flux_m + flux_p)

        flux_m = -dudz_m
        flux_p = -dudz_p
        flux_z = HALF*(flux_m + flux_p)


        !
        ! Compute boundary average flux
        !
        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        !
        ! Integrate flux
        !
        call worker%integrate_boundary(iu, integrand)


    end subroutine compute
    !**************************************************************************************************




end module LD_boundary_diffusive_operator
