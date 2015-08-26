!>  Boundary condition module
!!      - contains instantiations of all boundary condition for dynamically creating boundary conditions at run-time
!!      
!!
!!  Registering boundary conditions
!!      - To register a boundary condition:
!!          1st: Import it's definition for use in the current module
!!          2nd: Declare an instance of the boundary condition
!!          3rd: Extend 'create_bc' to include a selection criteria for the boundary condition
!!          4th: Under the selection criteria in 'create_bc', include a statement for dynamically allocating the boundary condition
!!
!--------------------------------------------------------
module mod_bc
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use atype_bc,       only: bc_t

    ! IMPORT BOUNDARY CONDITIONS
    use bc_linearadvection_extrapolate,   only: linearadvection_extrapolate_t
    implicit none


    ! Instantiate boundary conditions so they can be sourced by the factory
    type(linearadvection_extrapolate_t) :: LINEARADVECTION_EXTRAPOLATE




contains



    !> Boundary condition factory
    !!      - Allocate the appropriate boundary condition based on the incoming string specification
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      string  Character string used to select the appropriate boundary condition
    !!  @param[inout]   bc      Allocatable boundary condition
    !--------------------------------------------------------------
    subroutine create_bc(string,bc)
        character(*),                   intent(in)      :: string
        class(bc_t),    allocatable,    intent(inout)   :: bc

        integer(ik) :: ierr

        select case (string)


            ! Linear Advection - extrapolation boundary condition
            case ('extrapolate_la','extrapolation_la','Extrapolate_la','Extrapolation_la')
                allocate(bc, source=LINEARADVECTION_EXTRAPOLATE, stat=ierr)
                if (ierr /= 0) call AllocationError


            ! Euler - extrapolation boundary condition
            case ('extrapolate_euler','extrapolation_euler','Extrapolate_euler','Extrapolation_euler')
                call signal(FATAL,"create_bc: Euler extrapolation boundary condition is not yet implemented")



            ! DEFAULT - ERROR
            case default
                call signal(FATAL,"create_bc: Boundary condition string was not recognized. Check that the boundary condition is registered in create_bc and that spelling is correct.")            



        end select

    end subroutine




end module mod_bc
