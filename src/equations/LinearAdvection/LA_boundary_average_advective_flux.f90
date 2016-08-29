module LA_boundary_average_advective_flux
#include <messenger.h>
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: ZERO,ONE,TWO,HALF, ME, NEIGHBOR

    use type_boundary_flux,         only: boundary_flux_t
    use type_chidg_worker,          only: chidg_worker_t
    use type_properties,            only: properties_t
    use DNAD_D

    use LA_properties,              only: LA_properties_t
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
    type, extends(boundary_flux_t), public :: LA_boundary_average_advective_flux_t


    contains

        procedure   :: compute

    end type LA_boundary_average_advective_flux_t
    !********************************************************************************



contains

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
        class(LA_boundary_average_advective_flux_t),    intent(in)      :: self
        type(chidg_worker_t),                           intent(inout)   :: worker
        class(properties_t),                            intent(inout)   :: prop



        real(rk)                                :: cx, cy, cz
        integer(ik)                             :: iu


        type(AD_D), allocatable, dimension(:)   ::  &
            u_l, u_r, flux_x, flux_y, flux_z, integrand

        real(rk),   allocatable, dimension(:)   ::  &
            normx, normy, normz


        !
        ! Get variable index
        !
        iu = prop%get_eqn_index("u")



        !
        ! Get equation set properties
        !
        select type(prop)
            type is (LA_properties_t)
                cx = prop%c(1)
                cy = prop%c(2)
                cz = prop%c(3)
        end select

        
        !
        ! Interpolate solution to quadrature nodes
        !
        u_r = worker%interpolate(iu, 'value', ME)
        u_l = worker%interpolate(iu, 'value', NEIGHBOR)


        normx = worker%normal(1)
        normy = worker%normal(2)
        normz = worker%normal(3)


        !
        ! Compute boundary average flux
        !
        flux_x = ((cx*u_r + cx*u_l)/TWO ) * normx
        flux_y = ((cy*u_r + cy*u_l)/TWO ) * normy
        flux_z = ((cz*u_r + cz*u_l)/TWO ) * normz

        integrand = flux_x + flux_y + flux_z


        !
        ! Integrate flux
        !
        call worker%integrate_boundary(iu,integrand)


    end subroutine compute
    !**************************************************************************************************




end module LA_boundary_average_advective_flux
