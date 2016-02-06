module type_bcproperty
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: REQUIRED, OPTIONAL
    use type_function,  only: function_t
    use mod_function,   only: create_function
    implicit none
    private



    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !----------------------------------------------------------------------------------------
    type, public :: bcproperty_t

        character(len=:),   allocatable :: name_    !< Name of the property the function is applied to. Ex. Static Pressure.
        character(len=:),   allocatable :: type_    !< Property type.  'Required' or 'Optional'
        class(function_t),  allocatable :: fcn      !< allocatable function class

    contains

        procedure   :: set
        procedure   :: get
        procedure   :: status   !< Indicates if self%fcn component has been allocated

        procedure   :: get_name

        procedure   :: get_noptions
        procedure   :: get_option_key
        procedure   :: get_option_value

    end type bcproperty_t
    !***************************************************************************************




contains





    !>  Set boundary condition function information
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  comp    String indicating the component to set
    !!  @param[in]  str     String indicating the value to be set, or used to set.
    !!
    !---------------------------------------------------------------------------------------
    subroutine set(self,comp,str)
        class(bcproperty_t),   intent(inout)   :: self
        character(*),           intent(in)     :: comp
        character(*),           intent(in)     :: str



        select case (trim(comp))
            case ('Name','name')
                self%name_   = str

            case ('Type','type')
                self%type_   = str

            case ('Function','function','fcn')
                call create_function(self%fcn,str)

            case default
                call chidg_signal_one(FATAL,"bcfunction%set: Unrecognized component string.",comp)
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
        class(bcproperty_t),   intent(inout)   :: self
        character(*),           intent(in)      :: comp_str

        character(len=:), allocatable   :: comp_val


        select case (trim(comp_str))
            case ('Name','name')
                comp_val = self%name_


            case ('Type','type')
                comp_val = self%type_


            case ('Function','function')
                if ( allocated(self%fcn) ) then
                    comp_val = self%fcn%get_name()
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
        class(bcproperty_t),    intent(in)  :: self

        logical :: fcn_allocated
        

        fcn_allocated = allocated(self%fcn)
    

    end function status
    !****************************************************************************************







    
    !>  Return the name of the property.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------
    function get_name(self) result(name_string)
        class(bcproperty_t),    intent(in)  :: self

        character(len=:),   allocatable :: name_string

        name_string = self%name_

    end function get_name
    !***************************************************************************************













    !>  Return the number of available options.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_noptions(self) result(noptions)
        class(bcproperty_t),    intent(inout)  :: self
        
        integer(ik) :: noptions
        logical     :: fcn_allocated


        !
        ! Check that function is allocated
        !
        fcn_allocated = self%status()

        if ( fcn_allocated ) then

            noptions = self%fcn%get_noptions()

        else
            call chidg_signal(FATAL,"bcproperty%get_noptions: function not allocated")
        end if

    end function get_noptions
    !****************************************************************************************








    !>  Return an option key, given an option index.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_option_key(self,iopt) result(key)
        class(bcproperty_t),    intent(inout)  :: self
        integer(ik),            intent(in)  :: iopt

        character(len=:),   allocatable :: key
        logical     :: fcn_allocated

        !
        ! Check that function is allocated
        !
        fcn_allocated = self%status()

        if ( fcn_allocated ) then
            key = self%fcn%get_option_key(iopt)
        else
            call chidg_signal(FATAL,"bcproperty%get_option_key: function not allocated")
        end if

    end function get_option_key
    !****************************************************************************************








    !>  Return the value of an option, given an option key.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_option_value(self,key) result(val)
        class(bcproperty_t),    intent(inout)  :: self
        character(*),           intent(in)  :: key

        real(rk)        :: val
        logical     :: fcn_allocated

        !
        ! Check that function is allocated
        !
        fcn_allocated = self%status()

        if ( fcn_allocated ) then
            val = self%fcn%get_option_value(key)
        else
            call chidg_signal(FATAL,"bcproperty%get_option_value: function not allocated")
        end if


    end function get_option_value
    !****************************************************************************************



















end module type_bcproperty
