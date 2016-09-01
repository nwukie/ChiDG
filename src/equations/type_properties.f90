module type_properties
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use type_equation,  only: equation_t
    use type_fluid,     only: fluid_t
    implicit none




    !> Base properties type for storing equations, material definitions, 
    !! and miscelaneous data pertaining to a particular equationset
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------
    type, public :: properties_t
        
        ! Equations
        type(equation_t),   allocatable :: eqns(:)

        ! Materials
        class(fluid_t),     allocatable :: fluid

    contains

        procedure   :: nequations
        procedure   :: get_equation_index
        procedure   :: add_fluid

    end type properties_t
    !*********************************************************************************************






contains





    !> Search for a equation string in the self%eqns list. If found, return equation index.
    !! A set of equations could be stored in any order. So, when an equation is initialized, it
    !! is initialized with an index indicating its location in the set. That index is used to 
    !! access the correct solution data values.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2016
    !!
    !!  @param[in]  varstring   Character string identifying the desired variable
    !!
    !---------------------------------------------------------------------------------------------------
    function get_equation_index(self,varstring) result(varindex)
        class(properties_t),    intent(in)  :: self
        character(*),           intent(in)  :: varstring

        integer(ik) :: varindex, ieq
        logical     :: found = .false.

        varindex = 0


        ! Search for character string in self%eqns array. If found set index
        if (allocated(self%eqns)) then
            do ieq = 1,size(self%eqns)
                if (varstring == self%eqns(ieq)%name) then
                    varindex = self%eqns(ieq)%ind
                    found = .true.
                    exit
                end if
            end do

        else
            varindex = 0
        end if


!        ! Check if index was found
!        if (.not. found) call chidg_signal(FATAL,"Equation string not found in equation set properties")

    end function get_equation_index
    !***************************************************************************************************











    !>
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   8/31/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------------
    function nequations(self) result(neqns)
        class(properties_t),    intent(in)  :: self

        integer(ik) :: neqns

        if (allocated(self%eqns)) then
            neqns = size(self%eqns)
        else
            neqns = 0
        end if

    end function nequations
    !****************************************************************************************************















    !> Add a fluid definition to the properties type
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2016
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine add_fluid(self,fluid)
        class(properties_t),    intent(inout)   :: self
        class(fluid_t),         intent(in)      :: fluid

        integer(ik) :: ierr


        if (allocated(self%fluid)) deallocate(self%fluid)


        !
        ! Allocate new material definition
        !
        allocate(self%fluid, source=fluid, stat=ierr)
        if (ierr /= 0) call AllocationError

    end subroutine add_fluid
    !***************************************************************************************************







end module type_properties
