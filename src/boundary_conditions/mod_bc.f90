!>  Boundary condition module
!!      - contains instantiations of all boundary condition for dynamically creating boundary conditions at run-time
!!      
!!
!!  Registering boundary conditions
!!      - To register a boundary condition:
!!          1st: Import it's definition for use in the current module
!!          2nd: In the 'register_bcs' routine, declare an instance of the boundary condition
!!          3rd: In the 'register_bcs' routine, push an instanace of the bc to the registered_bcs vector 
!!
!--------------------------------------------------------
module mod_bc
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use type_bc_state,      only: bc_state_t
    use type_bcvector,      only: bcvector_t

    !
    ! Import boundary conditions
    !
    use bc_empty,                               only: empty_t
    use bc_periodic,                            only: periodic_t

    ! Scalar boundary conditions
    use bc_state_linearadvection_extrapolate,   only: linearadvection_extrapolate_t
    use bc_state_scalar_value,                  only: scalar_value_t
    use bc_state_scalar_extrapolate,            only: scalar_extrapolate_t


    ! Fluid boundary conditions
    use bc_state_wall,                          only: wall_t
    use bc_state_totalinlet,                    only: totalinlet_t
    use bc_state_pressureoutlet,                only: pressureoutlet_t



!    use bc_kirchoff,                        only: kirchoff_t

    implicit none


    !
    ! Global vector of registered boundary conditions
    !
    type(bcvector_t)    :: registered_bcs
    logical             :: initialized = .false.

contains


    !>  Register boundary conditions in a module vector.
    !!
    !!  This allows the available boundary conditions to be queried in the same way that they 
    !!  are registered for allocation. 
    !!
    !!  This gets called by chidg%init('env')
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine register_bcs()
        integer :: nbcs, ibc

        !
        ! Instantiate bcs
        !
        type(empty_t)                           :: EMPTY
        type(periodic_t)                        :: PERIODIC
        type(linearadvection_extrapolate_t)     :: LINEARADVECTION_EXTRAPOLATE
        type(scalar_value_t)                    :: SCALAR_VALUE
        type(scalar_extrapolate_t)              :: SCALAR_EXTRAPOLATE

        type(wall_t)                            :: WALL
        type(totalinlet_t)                      :: TOTALINLET
        type(pressureoutlet_t)                  :: PRESSUREOUTLET



        if ( .not. initialized ) then
            !
            ! Register in global vector
            !
            call registered_bcs%push_back(EMPTY)
            call registered_bcs%push_back(PERIODIC)

            call registered_bcs%push_back(LINEARADVECTION_EXTRAPOLATE)
            call registered_bcs%push_back(SCALAR_VALUE)
            call registered_bcs%push_back(SCALAR_EXTRAPOLATE)

            call registered_bcs%push_back(WALL)
            call registered_bcs%push_back(TOTALINLET)
            call registered_bcs%push_back(PRESSUREOUTLET)


            !
            ! Initialize each boundary condition in set. Doesn't need modified.
            !
            nbcs = registered_bcs%size()
            do ibc = 1,nbcs

                call registered_bcs%data(ibc)%state%init()

            end do

            !
            ! Confirm initialization
            !
            initialized = .true.

        end if


    end subroutine register_bcs
    !********************************************************************************************








    !>  Boundary condition factory
    !!      - Allocate a concrete boundary condition type based on the incoming string specification.
    !!      - Initialize the allocated boundary condition.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   1/31/2016
    !!
    !!  @param[in]      string  Character string used to select the appropriate boundary condition
    !!  @param[inout]   bc      Allocatable boundary condition
    !!
    !----------------------------------------------------------------------------------------------------
    subroutine create_bc(bcstring,bc)
        character(*),                       intent(in)      :: bcstring
        class(bc_state_t),  allocatable,    intent(inout)   :: bc

        integer(ik) :: ierr, bcindex


        if ( allocated(bc) ) then
            deallocate(bc)
        end if



        ! Find boundary condition string in 'registered_bcs' vector
        bcindex = registered_bcs%index_by_name(trim(bcstring))
        if (bcindex == 0) call chidg_signal_one(FATAL,"create_bc: boundary condition not recognized", trim(bcstring))


        ! Allocate conrete bc_t instance
        allocate(bc, source=registered_bcs%data(bcindex)%state, stat=ierr)
        if (ierr /= 0) call chidg_signal(FATAL,"create_bc: error allocating boundary condition from global vector.")


        ! Check boundary condition was allocated
        if ( .not. allocated(bc) ) call chidg_signal(FATAL,"create_bc: error allocating concrete boundary condition.")



    end subroutine create_bc
    !******************************************************************************************************







    !>  This is really just a utilitity for 'chidg edit' to dynamically list the avalable 
    !!  boundary conditions.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !-----------------------------------------------------------------------------------------------------
    subroutine list_bcs()
        integer                         :: nbcs, ibc
        character(len=:),   allocatable :: bcname

        nbcs = registered_bcs%size()


        do ibc = 1,nbcs

            bcname = registered_bcs%data(ibc)%state%get_name()
            call write_line(trim(bcname))

        end do ! ieqn

    end subroutine list_bcs
    !*****************************************************************************************************







    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   9/1/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------------------
    function check_bc_state_exists(state_string) result(state_found)
        character(len=*),   intent(in)  :: state_string

        integer(ik) :: state_index
        logical     :: state_found

        ! Find boundary condition string in 'registered_bcs' vector
        state_index = registered_bcs%index_by_name(trim(state_string))

        ! Set status of state_found
        state_found = (state_index /= 0)

    end function check_bc_state_exists
    !*******************************************************************************************************

































end module mod_bc
