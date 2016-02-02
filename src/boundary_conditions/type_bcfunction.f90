module type_bcfunction
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use type_function,  only: function_t
    implicit none
    private



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !----------------------------------------------------------------------------------------
    type, public :: bcfunction_t

        character(len=:),   allocatable :: name_    !< Name of the property the function is applied to. Ex. Static Pressure.
        character(len=:),   allocatable :: type_    !< Property type.  'Required' or 'Optional'

        class(function_t),  allocatable :: fcn      !< function

    contains

        procedure   :: set
        procedure   :: get
        procedure   :: status   !< Indicates if self%fcn component has been allocated

    end type bcfunction_t
    !***************************************************************************************




contains





    !>  Set boundary condition function information
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine set(self,comp_str,comp_val)
        class(bcfunction_t),   intent(inout)   :: self
        character(*),           intent(in)      :: comp_str
        character(*),           intent(in)      :: comp_val



        select case (trim(comp_str))
            case ('Name','name')
                self%name_   = comp_val

            case ('Type','type')
                self%type_   = comp_val

            case ('Function','function','fcn')
                call create_function(comp_val,self%fcn)

            case default
                call chidg_signal_one(FATAL,"bcfunction%set: Unrecognized component string.",comp_str)

        end select


    end subroutine set
    !****************************************************************************************











    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get(self,comp_str) result(comp_val)
        class(bcfunction_t),   intent(inout)   :: self
        character(*),           intent(in)      :: comp_str

        character(len=:), allocatable   :: comp_val


        select case (trim(comp_str))
            case ('Name','name')
                comp_val = self%name_



            case ('Type','type')
                comp_val = self%type_



            case ('Function','function')
                if ( allocated(self%fcn) ) then
                    comp_val = self%fcn%name
                else
                    call chidg_signal_one(WARN,"bcfunction%get: component not allocated.",comp_str)
                end if



            case default
                call chidg_signal_one(FATAL,"bcfunction%get: Unrecognized component string.",comp_str)

        end select



    end function get
    !***************************************************************************************








    !> Returns the status of the allocatable function component.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !----------------------------------------------------------------------------------------
    function status(self) result(fcn_allocated)
        class(bcfunction_t),    intent(in)  :: self

        logical :: fcn_allocated
        

        fcn_allocated = allocated(self%fcn)
    

    end function status
    !****************************************************************************************





end module type_bcfunction
