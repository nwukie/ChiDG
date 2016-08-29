module LD_boundary_diffusive_flux
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: ZERO,ONE,TWO,HALF, ME, NEIGHBOR

    use type_boundary_flux,         only: boundary_flux_t
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
    type, extends(boundary_flux_t), public :: LD_boundary_diffusive_flux_t


    contains

        procedure   :: compute

    end type LD_boundary_diffusive_flux_t
    !********************************************************************************



contains

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
        class(LD_boundary_diffusive_flux_t),    intent(in)      :: self
        type(chidg_worker_t),                   intent(inout)   :: worker
        class(properties_t),                    intent(inout)   :: prop


        integer(ik)                             :: iu

        type(AD_D), allocatable, dimension(:)   :: &
            dudx_m, dudx_p, flux_m, flux_p, flux, integrand

        real(rk),   allocatable, dimension(:)   :: &
            normx, normy, normz


        !
        ! Get variable index
        !
        iu = prop%get_eqn_index("u")



        !
        ! Interpolate solution to quadrature nodes
        !
        dudx_m = worker%interpolate(iu, 'ddx', ME)
        dudx_p = worker%interpolate(iu, 'ddx', NEIGHBOR)


        flux_m = -dudx_m
        flux_p = -dudx_p
        flux   = HALF*(flux_m + flux_p)

        !
        ! Compute boundary average flux
        !
        normx = worker%normal(1)
        integrand = flux * normx

        !
        ! Integrate flux
        !
        call worker%integrate_boundary(iu, integrand)


    end subroutine compute
    !**************************************************************************************************




end module LD_boundary_diffusive_flux
