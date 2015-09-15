module eqn_scalar
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ZERO,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use atype_equationset,      only: equationset_t
    use type_mesh,              only: mesh_t
    use atype_solverdata,       only: solverdata_t
    use type_properties,        only: properties_t
    use mod_interpolate,        only: interpolate
    use mod_integrate,          only: integrate_volume_flux, integrate_boundary_flux
    use mod_DNAD_tools,         only: compute_neighbor_face, compute_seed_element
    use DNAD_D

    use SCA_boundary_average_advective_flux, only: SCA_boundary_average_advective_flux_t
    use SCA_LaxFriedrichs_flux,              only: SCA_LaxFriedrichs_flux_t
    use SCA_volume_advective_flux,           only: SCA_volume_advective_flux_t
    use SCA_properties,                      only: SCA_properties_t
    implicit none

    private



    !   General scalar equation. Mostly used for unit/regression tests
    !   that need an initialized domain, so the scalar equation set can be used.
    !
    !   @author Nathan A. Wukie
    !
    !-----------------------------------------------------------
    type, extends(equationset_t), public :: scalar_e


    contains
        procedure   :: init

    end type scalar_e
    !-----------------------------------------------------------

contains


    !===========================================================
    !
    !   Equation set initialization
    !
    !===========================================================
    subroutine init(self)
        class(scalar_e), intent(inout) :: self


        type(SCA_volume_advective_flux_t)           :: volume_flux
        type(SCA_boundary_average_advective_flux_t) :: average_flux
        type(SCA_LaxFriedrichs_flux_t)              :: upwind_flux


        self%name = 'Scalar'

        !
        ! Allocate equation set properties
        !
        if (allocated(self%prop)) deallocate(self%prop)
        allocate(SCA_properties_t::self%prop)

        select type (prop => self%prop)
            type is (SCA_properties_t)
                prop%c(1) = ONE
                prop%c(2) = ZERO
                prop%c(3) = ZERO
        end select



        !
        ! Allocate and initialize equations
        !
        call self%add_equation("u",1)     


        !
        ! Add flux components
        !
        call self%add_volume_advective_flux(volume_flux)
        call self%add_boundary_advective_flux(average_flux)
        call self%add_boundary_advective_flux(upwind_flux)

    end subroutine







end module eqn_scalar
