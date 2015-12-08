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
    use type_bc,        only: bc_t

    ! IMPORT BOUNDARY CONDITIONS
    use bc_periodic,                    only: periodic_t
    use bc_linearadvection_extrapolate, only: linearadvection_extrapolate_t
    use bc_euler_wall,                  only: euler_wall_t
    use bc_euler_totalinlet,            only: euler_totalinlet_t
    use bc_euler_pressureoutlet,        only: euler_pressureoutlet_t
    use bc_euler_extrapolate,           only: euler_extrapolate_t
    use bc_lineuler_extrapolate,        only: lineuler_extrapolate_t
    use bc_lineuler_inlet,              only: lineuler_inlet_t
    implicit none


    ! Instantiate boundary conditions so they can be sourced by the factory
    type(periodic_t)                    :: PERIODIC
    type(linearadvection_extrapolate_t) :: LINEARADVECTION_EXTRAPOLATE
    type(euler_wall_t)                  :: EULER_WALL
    type(euler_totalinlet_t)            :: EULER_TOTALINLET
    type(euler_pressureoutlet_t)        :: EULER_PRESSUREOUTLET
    type(euler_extrapolate_t)           :: EULER_EXTRAPOLATE
    type(lineuler_extrapolate_t)        :: LINEULER_EXTRAPOLATE
    type(lineuler_inlet_t)              :: LINEULER_INLET



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

        select case (trim(string))
            ! PERIODIC
            case ('periodic','Periodic')
                allocate(bc, source=PERIODIC, stat=ierr)

            ! Linear Advection - extrapolation boundary condition
            case ('extrapolate_la','extrapolation_la','Extrapolate_la','Extrapolation_la')
                allocate(bc, source=LINEARADVECTION_EXTRAPOLATE, stat=ierr)


            ! Euler - extrapolation boundary condition
            case ('extrapolate_euler','extrapolation_euler','Extrapolate_euler','Extrapolation_euler')
                call chidg_signal(FATAL,"create_bc: Euler extrapolation boundary condition is not yet implemented")

            ! Euler - slip wall
            case ('euler_wall','slip_wall','Euler_Wall','Slip_Wall','SlipWall')
                allocate(bc, source=EULER_WALL, stat=ierr)

            ! Euler - total inlet
            case ('euler_totalinlet','Euler_TotalInlet')
                allocate(bc, source=EULER_TOTALINLET, stat=ierr)



            ! Euler - pressure oulet
            case ('euler_pressureoutlet','Euler_PressureOutlet')
                allocate(bc, source=EULER_PRESSUREOUTLET, stat=ierr)



            ! Euler - extrapolate
            case ('euler_extrapolate','Euler_Extrapolate')
                allocate(bc, source=EULER_EXTRAPOLATE, stat=ierr)



            ! Linearized Euler
            case ('lineuler_extrapolate')
                allocate(bc, source=LINEULER_EXTRAPOLATE, stat=ierr)


            case ('lineuler_inlet')
                allocate(bc, source=LINEULER_INLET, stat=ierr)





            ! DEFAULT - ERROR
            case default
                call chidg_signal(FATAL,"create_bc: Boundary condition string was not recognized. Check that the boundary condition is registered in create_bc and that spelling is correct.")            



        end select


        !
        ! Check allocation status
        !
        if (ierr /= 0) call AllocationError


    end subroutine




end module mod_bc
