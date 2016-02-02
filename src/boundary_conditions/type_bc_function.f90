module type_bc_function
#include <messenger.h>
    implicit none
    private



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !----------------------------------------------------------------------------------------
    type, public :: bc_function_t

        character(len=:),   allocatable :: name_    !< Name of the property the function is applied to. Ex. Static Pressure.
        character(len=:),   allocatable :: type_    !< Property type.  'Required' or 'Optional'
        character(len=:),   allocatable :: status_  !< Property status. 'Set', 'Empty'

        class(function_t),  allocatable :: fcn

    contains

        procedure   :: set
        procedure   :: get

    end type bc_function_t
    !***************************************************************************************




contains





    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine set(self,comp_str,comp_val)
        class(bc_function_t),   intent(inout)   :: self
        character(*),           intent(in)      :: comp_str
        character(*),           intent(in)      :: comp_val



        select case (trim(comp_str))
            case ('Name','name')
                self%name_   = comp_val

            case ('Type','type')
                self%type_   = comp_val

            case ('Status','status')
                self%status_ = comp_val

            case default
                call chidg_signal_one(FATAL,"bc_function%set: Unrecognized component string.",comp_str)

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
        class(bc_function_t),   intent(inout)   :: self
        character(*),           intent(in)      :: comp_str

        character(len=:), allocatable   :: comp_val


        select case (trim(comp_str))
            case ('Name','name')
                comp_val = self%name_

            case ('Type','type')
                comp_val = self%type_

            case ('Status','status')
                comp_val = self%status_

            case default
                call chidg_signal_one(FATAL,"bc_function%set: Unrecognized component string.",comp_str)

        end select



    end function get
    !***************************************************************************************











end module type_bc_function
