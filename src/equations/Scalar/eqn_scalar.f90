module eqn_scalar
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,ONE,TWO

    use type_equationset,       only: equationset_t
    use type_properties,        only: properties_t

    use SCA_boundary_average_advective_flux, only: SCA_boundary_average_advective_flux_t
    use SCA_LaxFriedrichs_flux,              only: SCA_LaxFriedrichs_flux_t
    use SCA_volume_advective_flux,           only: SCA_volume_advective_flux_t
    use SCA_properties,                      only: SCA_properties_t
    implicit none

    private





    !>  General scalar equation. Mostly used for unit/regression tests
    !!  that need an initialized domain, so the scalar equation set can be used.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    type, extends(equationset_t), public :: scalar_e


    contains

        procedure   :: init

    end type scalar_e
    !****************************************************************************************************





contains





    !> Initialize scalar equation set.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine init(self)
        class(scalar_e), intent(inout) :: self

        type(SCA_volume_advective_flux_t)           :: volume_flux
        type(SCA_boundary_average_advective_flux_t) :: average_flux
        type(SCA_LaxFriedrichs_flux_t)              :: upwind_flux
        type(SCA_properties_t)                      :: prop


        !
        ! Set equationset name.
        !
        call self%set_name("Scalar")


        !
        ! Set up properties.
        !
        prop%c(1) = ONE
        prop%c(2) = ZERO
        prop%c(3) = ZERO


        !
        ! Allocate equation set properties.
        !
        call self%add_properties(prop)


        !
        ! Allocate and initialize equations.
        !
        call self%add_equation("u",1)     


        !
        ! Add flux components.
        !
        call self%add_volume_advective_flux(volume_flux)
        call self%add_boundary_advective_flux(average_flux)
        call self%add_boundary_advective_flux(upwind_flux)


    end subroutine init
    !****************************************************************************************************







end module eqn_scalar
