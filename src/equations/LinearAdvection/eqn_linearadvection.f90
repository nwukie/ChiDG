module eqn_linearadvection
#include <messenger.h>
    use mod_kinds,                          only: rk,ik
    use mod_constants,                      only: NFACES,ZERO,ONE,TWO,HALF, &
                                                  XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use type_equationset,                   only: equationset_t
    use type_properties,                    only: properties_t

    use LA_boundary_average_advective_flux, only: LA_boundary_average_advective_flux_t
    use LA_LaxFriedrichs_flux,              only: LA_LaxFriedrichs_flux_t
    use LA_volume_advective_flux,           only: LA_volume_advective_flux_t
    use LA_properties,                      only: LA_properties_t
    implicit none

    private



    !>
    !!
    !!
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------------
    type, extends(equationset_t), public :: linearadvection_e


    contains

        procedure   :: init

    end type linearadvection_e
    !******************************************************************************************************






contains




    !>
    !!
    !!
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------------
    subroutine init(self)
        class(linearadvection_e), intent(inout) :: self


        type(LA_boundary_average_advective_flux_t)  :: average_flux
        type(LA_LaxFriedrichs_flux_t)               :: upwind_flux
        type(LA_volume_advective_flux_t)            :: volume_flux
        type(LA_properties_t)                       :: prop

        
        !
        ! Set equationset name.
        !
        call self%set_name("LinearAdvection")


        !
        ! Set properties
        !
        prop%c(1) = ONE
        prop%c(2) = ZERO
        prop%c(3) = ZERO


        !
        ! Allocate equation set properties
        !
        call self%add_properties(prop)


        !
        ! Allocate and initialize equations
        !
        call self%add_equation("u",1)     



        !
        ! Allocate flux components to specific types for the equation set
        ! 
        call self%add_volume_advective_flux(volume_flux)
        call self%add_boundary_advective_flux(average_flux)
        call self%add_boundary_advective_flux(upwind_flux)



    end subroutine init
    !*********************************************************************************************************





end module eqn_linearadvection
