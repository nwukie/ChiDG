module SD_boundary_operator
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: ZERO,ONE,TWO,HALF
    use type_operator,              only: operator_t
    use type_chidg_worker,          only: chidg_worker_t
    use type_properties,            only: properties_t
    use DNAD_D

    use SD_properties,              only: SD_properties_t
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
    type, extends(operator_t), public :: SD_boundary_operator_t


    contains

        procedure   :: init
        procedure   :: compute

    end type SD_boundary_operator_t
    !********************************************************************************



contains


    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(SD_boundary_operator_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Scalar Diffusion Boundary Average Operator")

        !
        ! Set operator type
        !
        call self%set_operator_type("Boundary Diffusive Operator")

        !
        ! Set operator equations
        !
        call self%add_primary_field("u")

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
        class(SD_boundary_operator_t),  intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop



        type(AD_D), allocatable, dimension(:)   :: &
            dudx_m, dudy_m, dudz_m, dudx_p, dudy_p, dudz_p, &
            flux_x, flux_y, flux_z, flux_m, flux_p, flux, integrand, &
            mu_m, mu_p

        real(rk),   allocatable, dimension(:)   :: &
            normx, normy, normz


        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)


        !
        ! Interpolate solution to quadrature nodes
        !
        dudx_m = worker%get_primary_field_face('u', 'ddx + lift', 'face interior')
        dudy_m = worker%get_primary_field_face('u', 'ddy + lift', 'face interior')
        dudz_m = worker%get_primary_field_face('u', 'ddz + lift', 'face interior')


        dudx_p = worker%get_primary_field_face('u', 'ddx + lift', 'face exterior')
        dudy_p = worker%get_primary_field_face('u', 'ddy + lift', 'face exterior')
        dudz_p = worker%get_primary_field_face('u', 'ddz + lift', 'face exterior')


        !
        ! Compute scalar coefficient
        !
        mu_m = worker%get_model_field_face('Scalar Diffusion Coefficient', 'value', 'face interior')
        mu_p = worker%get_model_field_face('Scalar Diffusion Coefficient', 'value', 'face exterior')


        flux_m = -mu_m*dudx_m
        flux_p = -mu_p*dudx_p
        flux_x = HALF*(flux_m + flux_p)

        flux_m = -mu_m*dudy_m
        flux_p = -mu_p*dudy_p
        flux_y = HALF*(flux_m + flux_p)

        flux_m = -mu_m*dudz_m
        flux_p = -mu_p*dudz_p
        flux_z = HALF*(flux_m + flux_p)


        !
        ! Compute boundary average flux
        !
        integrand = flux_x*normx + flux_y*normy + flux_z*normz

        !
        ! Integrate flux
        !
        call worker%integrate_boundary('u',integrand)


    end subroutine compute
    !**************************************************************************************************




end module SD_boundary_operator
