module eqn_duallinearadvection
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NFACES,ZERO,ONE,TWO,HALF, &
                                      XI_MIN,XI_MAX,ETA_MIN,ETA_MAX,ZETA_MIN,ZETA_MAX,DIAG

    use type_equationset,       only: equationset_t
    use type_mesh,              only: mesh_t
    use type_solverdata,        only: solverdata_t

    use DLA_boundary_average_advective_flux,    only: DLA_boundary_average_advective_flux_t
    use DLA_LaxFriedrichs_flux,                 only: DLA_LaxFriedrichs_flux_t
    use DLA_volume_advective_flux,              only: DLA_volume_advective_flux_t
    use DLA_properties,                         only: DLA_properties_t
    implicit none

    private





    !> This equation set exists really just to test equationsets with 
    !! more than one equation. The idea is just to compute the linear
    !! advecdtion solution twice at the same time. The equations are 
    !! independent of each other. So, we can verify, for example,
    !! the volume flux jacobians for each equation. They should be the
    !! same as for the single LinearAdvection equation set
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    type, extends(equationset_t), public :: duallinearadvection_e


    contains

        procedure   :: init

    end type duallinearadvection_e
    !*************************************************************************************************





contains


    !>
    !!
    !!
    !!
    !!
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine init(self)
        class(duallinearadvection_e), intent(inout) :: self

        type(DLA_volume_advective_flux_t)           :: volume_flux
        type(DLA_boundary_average_advective_flux_t) :: average_flux
        type(DLA_LaxFriedrichs_flux_t)              :: upwind_flux
        type(DLA_properties_t)                      :: prop


        !
        ! Set equationset name.
        !
        self%name = 'DualLinearAdvection'


        !
        ! Allocate equation set properties
        !
        call self%add_properties(prop)


        !
        ! Allocate and initialize equations
        !
        call self%add_equation('u_a',1)
        call self%add_equation('u_b',2)



        !
        ! Attach flux components
        !
        call self%add_volume_advective_flux(volume_flux)
        call self%add_boundary_advective_flux(average_flux)
        call self%add_boundary_advective_flux(upwind_flux)


    end subroutine init
    !*************************************************************************************************








end module eqn_duallinearadvection
