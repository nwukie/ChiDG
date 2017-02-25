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
    use bc_state_scalar_derivative,             only: scalar_derivative_t
    use bc_state_scalar_extrapolate,            only: scalar_extrapolate_t


    ! Fluid boundary conditions
    use bc_state_wall,                          only: wall_t
    use bc_state_moving_wall,                   only: moving_wall_t
    use bc_state_inlet_total,                   only: inlet_total_t
    use bc_state_outlet_constant_pressure,      only: outlet_constant_pressure_t
    use bc_state_outlet_point_pressure,         only: outlet_point_pressure_t
    use bc_state_fluid_extrapolate,             only: fluid_extrapolate_t
    use bc_state_momentum_inlet,                only: momentum_inlet_t
    use bc_state_symmetry,                      only: symmetry_t
    use bc_state_farfield,                      only: farfield_t

    ! Turbulence boundary conditions
    use bc_state_spalart_allmaras_inlet,        only: spalart_allmaras_inlet_t
    use bc_state_spalart_allmaras_outlet,       only: spalart_allmaras_outlet_t
    use bc_state_spalart_allmaras_symmetry,     only: spalart_allmaras_symmetry_t
    use bc_state_spalart_allmaras_farfield,     only: spalart_allmaras_farfield_t
    use bc_state_spalart_allmaras_wall,         only: spalart_allmaras_wall_t

    ! Artificial Viscosity boundary conditions
    use bc_state_artificial_viscosity_wall,     only: artificial_viscosity_wall_t
    use bc_state_artificial_viscosity_inlet,    only: artificial_viscosity_inlet_t
    use bc_state_artificial_viscosity_outlet,   only: artificial_viscosity_outlet_t
    use bc_state_artificial_viscosity_symmetry, only: artificial_viscosity_symmetry_t


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
        type(scalar_derivative_t)               :: SCALAR_DERIVATIVE
        type(scalar_extrapolate_t)              :: SCALAR_EXTRAPOLATE

        type(wall_t)                            :: WALL
        type(moving_wall_t)                     :: MOVING_WALL
        type(inlet_total_t)                     :: INLET_TOTAL
        type(outlet_constant_pressure_t)        :: OUTLET_CONSTANT_PRESSURE
        type(outlet_point_pressure_t)           :: OUTLET_POINT_PRESSURE
        type(fluid_extrapolate_t)               :: FLUID_EXTRAPOLATE
        type(momentum_inlet_t)                  :: MOMENTUM_INLET
        type(symmetry_t)                        :: SYMMETRY
        type(farfield_t)                        :: FARFIELD

        type(spalart_allmaras_inlet_t)          :: SPALART_ALLMARAS_INLET
        type(spalart_allmaras_outlet_t)         :: SPALART_ALLMARAS_OUTLET
        type(spalart_allmaras_symmetry_t)       :: SPALART_ALLMARAS_SYMMETRY
        type(spalart_allmaras_farfield_t)       :: SPALART_ALLMARAS_FARFIELD
        type(spalart_allmaras_wall_t)           :: SPALART_ALLMARAS_WALL

        type(artificial_viscosity_wall_t)       :: ARTIFICIAL_VISCOSITY_WALL
        type(artificial_viscosity_inlet_t)      :: ARTIFICIAL_VISCOSITY_INLET
        type(artificial_viscosity_outlet_t)     :: ARTIFICIAL_VISCOSITY_OUTLET
        type(artificial_viscosity_symmetry_t)   :: ARTIFICIAL_VISCOSITY_SYMMETRY



        if ( .not. initialized ) then
            !
            ! Register in global vector
            !
            call registered_bcs%push_back(EMPTY)
            call registered_bcs%push_back(PERIODIC)

            call registered_bcs%push_back(LINEARADVECTION_EXTRAPOLATE)
            call registered_bcs%push_back(SCALAR_VALUE)
            call registered_bcs%push_back(SCALAR_DERIVATIVE)
            call registered_bcs%push_back(SCALAR_EXTRAPOLATE)

            call registered_bcs%push_back(WALL)
            call registered_bcs%push_back(MOVING_WALL)
            call registered_bcs%push_back(INLET_TOTAL)
            call registered_bcs%push_back(OUTLET_CONSTANT_PRESSURE)
            call registered_bcs%push_back(OUTLET_POINT_PRESSURE)
            call registered_bcs%push_back(FLUID_EXTRAPOLATE)
            call registered_bcs%push_back(MOMENTUM_INLET)
            call registered_bcs%push_back(SYMMETRY)
            call registered_bcs%push_back(FARFIELD)

            call registered_bcs%push_back(SPALART_ALLMARAS_INLET)
            call registered_bcs%push_back(SPALART_ALLMARAS_OUTLET)
            call registered_bcs%push_back(SPALART_ALLMARAS_SYMMETRY)
            call registered_bcs%push_back(SPALART_ALLMARAS_FARFIELD)
            call registered_bcs%push_back(SPALART_ALLMARAS_WALL)

            call registered_bcs%push_back(ARTIFICIAL_VISCOSITY_WALL)
            call registered_bcs%push_back(ARTIFICIAL_VISCOSITY_INLET)
            call registered_bcs%push_back(ARTIFICIAL_VISCOSITY_OUTLET)
            call registered_bcs%push_back(ARTIFICIAL_VISCOSITY_SYMMETRY)



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
        integer                     :: nbcs, ibc
        character(:),   allocatable :: bcname

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
    function check_bc_state_registered(state_string) result(state_found)
        character(len=*),   intent(in)  :: state_string

        integer(ik) :: state_index
        logical     :: state_found

        ! Find boundary condition string in 'registered_bcs' vector
        state_index = registered_bcs%index_by_name(trim(state_string))

        ! Set status of state_found
        state_found = (state_index /= 0)

    end function check_bc_state_registered
    !*******************************************************************************************************

































end module mod_bc
