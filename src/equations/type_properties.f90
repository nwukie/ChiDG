module type_properties
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use type_equation,  only: equation_t
    use atype_fluid,    only: fluid_t
    use atype_solid,    only: solid_t

    implicit none

    private



    !> Base properties type for storing equations, material definitions, 
    !! and miscelaneous data pertaining to a particular equationset
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------
    type, public :: properties_t
        
        ! Equations
        type(equation_t), allocatable   :: eqns(:)


        ! Materials
        class(fluid_t),   allocatable   :: fluid
        class(solid_t),   allocatable   :: solid

    contains
        procedure   :: init
        procedure   :: get_eqn_index

    end type properties_t
    !-----------------------------------------------------------------






contains

    !>  Equation properties initialization
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  fluid   fluid_t object to be assigned
    !!  @param[in]  solid   solid_t object to be assigned
    !!
    !------------------------------------------------------------------
    subroutine init(self,fluid,solid)
        class(properties_t),    intent(inout)            :: self
        class(fluid_t),         intent(in),    optional  :: fluid
        class(solid_t),         intent(in),    optional  :: solid

        integer(ik) :: ierr

        if (present(fluid)) then
            allocate(self%fluid, source=fluid, stat=ierr)
            if (ierr /= 0) call AllocationError
        end if

        if (present(solid)) then
            allocate(self%solid, source=solid, stat=ierr)
            if (ierr /= 0) call AllocationError
        end if

    end subroutine








    !> Search for a equation string in the self%eqns list. If found, return equation index.
    !! A set of equations could be stored in any order. So, when an equation is initialized, it
    !! is initialized with an index indicating its location in the set. That index is used to 
    !! access the correct solution data values.
    !!
    !!   @author Nathan A. Wukie
    !!
    !!   @param[in]  varstring   Character string identifying the desired variable
    !!
    !-------------------------------------------------------------------------------------
    function get_eqn_index(self,varstring) result(varindex)
        class(properties_t),    intent(in)  :: self
        character(*),           intent(in)  :: varstring

        integer(ik) :: varindex, ieq
        logical     :: found = .false.

        varindex = 123456789

        !
        ! Search for character string in self%eqns array. If found set index
        !
        do ieq = 1,size(self%eqns)
            if (varstring == self%eqns(ieq)%name) then
                varindex = self%eqns(ieq)%ind
                found = .true.
                exit
            end if
        end do

        !
        ! Check if index was found
        !
        if (.not. found) call chidg_signal(FATAL,"Equation string not found in equation set properties")

    end function









end module type_properties
