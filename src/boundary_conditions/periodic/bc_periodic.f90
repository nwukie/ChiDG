module bc_periodic
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ONE, TWO, HALF, LOCAL
    use type_bc,            only: bc_t
    use type_solverdata,    only: solverdata_t
    use type_mesh,          only: mesh_t
    use type_properties,    only: properties_t
    use type_face_info,     only: face_info_t
    use type_function_info, only: function_info_t
    implicit none





    !> Periodic boundary condition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !--------------------------------------------------------------------------------
    type, extends(bc_t) :: periodic_t

    contains

        procedure add_options
        procedure compute

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
    subroutine compute(self,mesh,sdata,prop,face,flux)
        class(periodic_t),              intent(inout)   :: self
        type(mesh_t),                   intent(in)      :: mesh(:)
        type(solverdata_t),             intent(inout)   :: sdata
        class(properties_t),            intent(inout)   :: prop
        type(face_info_t),              intent(in)      :: face
        type(function_info_t),          intent(in)      :: flux


        ! DO NOTHING IN PERIODIC BOUNDARY CONDITION

    end subroutine compute
    !***********************************************************************************













end module bc_periodic
