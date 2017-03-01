module type_bc_group
#include <messenger.h>
    use type_bcvector,  only: bcvector_t
    use type_bc_state,  only: bc_state_t
    implicit none




    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/9/2016
    !!
    !-------------------------------------------------------------------------
    type, public :: bc_group_t
        
        character(:),       allocatable :: name         ! Boundary State Group name
        character(:),       allocatable :: family       ! Boundary State Group family
        type(bcvector_t)                :: bc_states    ! Vector of boundary condition state functions for each group.

    contains

        procedure   :: set_name
        procedure   :: get_name
        procedure   :: set_family
        procedure   :: get_family

        procedure   :: add_bc_state

    end type bc_group_t
    !*************************************************************************



contains



    !>  Set the bc_group name.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/1/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine set_name(self,name)
        class(bc_group_t),  intent(inout)   :: self
        character(*),       intent(in)      :: name

        self%name = name

    end subroutine set_name
    !******************************************************************************************



    !>  Return bc_group name.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/1/2017
    !!
    !------------------------------------------------------------------------------------------
    function get_name(self) result(name)
        class(bc_group_t),  intent(inout)   :: self

        character(:),   allocatable :: name, user_msg

        if (allocated(self%name)) then
            name = self%name
        else
            user_msg = "bc%get_name: It looks like the boundary condition group name was never set.&
                        Make sure bc_group%set_name gets called in the boundary condition group &
                        initialization routine"
            call chidg_signal(FATAL,user_msg)
        end if

    end function get_name
    !******************************************************************************************





    !>  Set the bc_group family.
    !!
    !!  bc_family may be:
    !!      - Wall
    !!      - Inlet
    !!      - Outlet
    !!      - Symmetry
    !!      - Periodic
    !!      - Farfield
    !!      - Scalar
    !!      - Extrapolation
    !!      - Empty
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/1/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine set_family(self,family)
        class(bc_group_t),  intent(inout)   :: self
        character(*),       intent(in)      :: family

        character(:),   allocatable :: user_msg

        if ( (trim(family) == 'Wall'    )       .or. &
             (trim(family) == 'Inlet'   )       .or. &
             (trim(family) == 'Outlet'  )       .or. &
             (trim(family) == 'Symmetry')       .or. &
             (trim(family) == 'Periodic')       .or. &
             (trim(family) == 'Farfield')       .or. &
             (trim(family) == 'Scalar'  )       .or. &
             (trim(family) == 'Extrapolation')  .or. &
             (trim(family) == 'Empty'   ) ) then

            self%family = family

        else
            user_msg = "bc_group%set_family: The string passed in to set the boundary condition family did &
                        not match any of valid boundary condition families. These include: 'Wall', &
                        'Inlet', 'Outlet', 'Symmetry', 'Periodic', 'Farfield', 'Scalar', 'Extrapolation'"
            call chidg_signal_one(FATAL,user_msg,family)
        end if

    end subroutine set_family
    !******************************************************************************************






    !>  Return the bc_group family.
    !!
    !!  bc_family may be:
    !!      'Wall', 'Inlet', 'Outlet', 'Symmetry', 'Periodic', 'Farfield', 'Scalar'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/1/2017
    !!
    !!
    !------------------------------------------------------------------------------------------
    function get_family(self) result(family)
        class(bc_group_t),  intent(inout)   :: self

        character(:),   allocatable :: family, user_msg

        if (allocated(self%family)) then
            family = self%family
        else
            family = ' '
        end if

    end function get_family
    !******************************************************************************************











    !>  Add a bc_state object to the bc_group.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/1/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine add_bc_state(self,bc_state)
        class(bc_group_t),  intent(inout)   :: self
        class(bc_state_t),  intent(in)      :: bc_state 

        character(:),   allocatable :: group_family, state_family, user_msg
        logical                     :: add_state
        

        !
        ! Check bc_state family conforms to any already added.
        !
        group_family = self%get_family()
        state_family = bc_state%get_family()
        add_state = (group_family == ' ') .or. (trim(group_family) == trim(state_family))


        !
        ! Add to vector of bc_states on the group.
        !
        if (add_state) then
            call self%set_family(trim(state_family))
            call self%bc_states%push_back(bc_state)
        else
            user_msg = "bc_group%add_bc_state: An attempt was made to add a bc_state object &
                        to a bc_group with dissimilar family. As a rule, bc_group objects &
                        may only contain bc_state objects of a single family."
            call chidg_signal_one(FATAL,user_msg, bc_state%get_name())
        end if

    end subroutine add_bc_state
    !******************************************************************************************









end module type_bc_group
