module bc_periodic
#include <messenger.h>
    use type_bc_state,      only: bc_state_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    implicit none





    !> Periodic boundary condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !--------------------------------------------------------------------------------
    type, extends(bc_state_t) :: periodic_t

    contains

        procedure   :: init
        procedure   :: compute_bc_state

    end type periodic_t
    !********************************************************************************



contains



    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/29/2016
    !!
    !--------------------------------------------------------------------------------
    subroutine init(self)
        class(periodic_t),   intent(inout) :: self
        
        !
        ! Set operator name
        !
        call self%set_name("Periodic")




        !
        ! Add functions
        !
        call self%bcproperties%add('type',          'Required')
        call self%bcproperties%add('offset_x',      'Required')
        call self%bcproperties%add('offset_y',      'Required')
        call self%bcproperties%add('offset_z',      'Required')
        call self%bcproperties%add('offset_theta',  'Required')


        !
        ! Add parameters
        !



    end subroutine init
    !********************************************************************************






    !> Boundary condition compute routine called by spatial scheme
    !!      - Matching periodic boundary condition, so the interior scheme 
    !!        is just adjusted to connect elements and no extra calculation
    !!        routine is needed here. Hence the empty routine below.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !-----------------------------------------------------------------------------------
    subroutine compute_bc_state(self,worker,prop)
        class(periodic_t),              intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop


        ! DO NOTHING IN PERIODIC BOUNDARY CONDITION

    end subroutine compute_bc_state
    !***********************************************************************************













end module bc_periodic
