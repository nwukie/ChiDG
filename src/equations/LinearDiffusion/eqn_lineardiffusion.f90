module eqn_lineardiffusion
#include <messenger.h>
    use mod_kinds,                      only: rk,ik
    use mod_constants,                  only: NFACES,ZERO,ONE,TWO,HALF, &
                                              XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use type_equationset,               only: equationset_t
    use type_properties,                only: properties_t

    use LD_boundary_diffusive_flux,     only: LD_boundary_diffusive_flux_t
    use LD_volume_diffusive_flux,       only: LD_volume_diffusive_flux_t
    use LD_properties,                  only: LD_properties_t
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
    type, extends(equationset_t), public :: lineardiffusion_e


    contains

        procedure   :: init

    end type lineardiffusion_e
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
        class(lineardiffusion_e), intent(inout) :: self


        type(LD_boundary_diffusive_flux_t)  :: boundary_flux
        type(LD_volume_diffusive_flux_t)    :: volume_flux
        type(LD_properties_t)               :: prop

        
        !
        ! Set equationset name.
        !
        call self%set_name("LinearDiffusion")


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
        call self%add_volume_diffusive_flux(volume_flux)
        call self%add_boundary_diffusive_flux(boundary_flux)



    end subroutine init
    !*********************************************************************************************************





end module eqn_lineardiffusion
