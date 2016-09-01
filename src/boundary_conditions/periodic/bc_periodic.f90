module bc_periodic
#include <messenger.h>
    use type_bc_operator,   only: bc_operator_t
    use type_chidg_worker,  only: chidg_worker_t
    use type_properties,    only: properties_t
    implicit none





    !> Periodic boundary condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !--------------------------------------------------------------------------------
    type, extends(bc_operator_t) :: periodic_t

    contains

        procedure   :: add_options
        procedure   :: init
        procedure   :: compute

    end type periodic_t
    !********************************************************************************



contains



    !>  Procedure for registering boundary condition options. Needs executed upon allocation.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_options(self)    
        class(periodic_t),  intent(inout)   :: self

        !
        ! Set name
        !
        call self%set_name('periodic')


        !
        ! Add functions
        !
        call self%bcproperties%add('type', 'Required')
        call self%bcproperties%add('offset_x', 'Required')
        call self%bcproperties%add('offset_y', 'Required')
        call self%bcproperties%add('offset_z', 'Required')
        call self%bcproperties%add('offset_theta', 'Required')


        !
        ! Add parameters
        !


    end subroutine add_options
    !******************************************************************************************




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
        call self%set_name("periodic")

        !
        ! Set operator type
        !
        call self%set_operator_type("Boundary Advective Flux")

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
    !subroutine compute(self,mesh,sdata,prop,face,fcn)
    subroutine compute(self,worker,prop)
        class(periodic_t),              intent(inout)   :: self
        type(chidg_worker_t),           intent(inout)   :: worker
        class(properties_t),            intent(inout)   :: prop
!        type(mesh_t),                   intent(in)      :: mesh(:)
!        type(solverdata_t),             intent(inout)   :: sdata
!        type(face_info_t),              intent(in)      :: face
!        type(function_info_t),          intent(in)      :: fcn


        ! DO NOTHING IN PERIODIC BOUNDARY CONDITION

    end subroutine compute
    !***********************************************************************************













end module bc_periodic
