module type_bcproperty
#include <messenger.h>
    use mod_kinds,      only: rk, ik
    use mod_constants,  only: REQUIRED, OPTIONAL
    use type_function,  only: function_t
    use mod_function,   only: create_function
    implicit none
    private



    !>  A class for boundary condition properties.
    !!
    !!  If one wishes to allow a user to specify information about a boundary condition, that
    !!  occurs via this boundary condition property mechanism. When implementing a boundary
    !!  condition, a developer decides what information should be specified at runtime by the user.
    !!  For example, total pressure and total temperature. The developer then adds an instance of
    !!  this bcproperty_t for each property. One for total pressure. One for total temperature.
    !!
    !!  The bcproperty_t class then manages the options for a given property. A function must be 
    !!  specified to describe the property in space and time. This could be a constant, or it
    !!  could vary in space and time, depending on the allocation of self%fcn. The dynamic
    !!  function_t then has potential options to be set. For a constant function, just a value.
    !!  For a trig function, maybe amplitude and phase. These particular options can be accessed
    !!  through the bcproperty_t methods get_noptions, get_option_key, get_option_value.
    !!  
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !----------------------------------------------------------------------------------------
    type, public :: bcproperty_t

        character(len=:),   allocatable :: name_    !< Name of the property the function is applied to. Ex. Static Pressure.
        character(len=:),   allocatable :: type_    !< Property type.  'Function' or 'Parameter'
        class(function_t),  allocatable :: fcn      !< allocatable function class

    contains

        procedure   :: status           !< Indicates if self%fcn component has been allocated

        procedure   :: set              !<  Set the name_, type_, or fcn components
        procedure   :: get              !<  Get the name_, type_, or fcn compoenents

        procedure   :: get_name         !< Return the name of the bcproperty. Ex. StaticPressure.
        procedure   :: get_noptions     !< Return the number of available options for the property. Depends on fcn.
        procedure   :: get_option_key   !< Return the key of an option, given option index.
        procedure   :: get_option_value !< Return the value of an option, given a key.

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
        class(bcproperty_t),   intent(in)   :: self
        character(*),          intent(in)   :: comp_str

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








    !> Returns the status of the allocatable function component. True if allocated, false otherwise.
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
    !---------------------------------------------------------------------------------------
    function get_name(self) result(name_string)
        class(bcproperty_t),    intent(in)  :: self

        character(len=:),   allocatable :: name_string

        name_string = self%name_

    end function get_name
    !***************************************************************************************













    !>  Return the number of available options. Depends on the allocated self%fcn
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function get_noptions(self) result(noptions)
        class(bcproperty_t),    intent(in)  :: self
        
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
    !!  @param[in]  iopt    Integer index of the option to be queried.
    !!  @result     key     String containing the name of the queried key.
    !!
    !----------------------------------------------------------------------------------------
    function get_option_key(self,iopt) result(key)
        class(bcproperty_t),    intent(in)  :: self
        integer(ik),            intent(in)  :: iopt

        character(len=:),   allocatable :: key
        logical                         :: fcn_allocated

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
    !!  @param[in]  key     String(key) indicating the option to be queried.
    !!  @result     val     Value of the specified key.
    !!
    !-----------------------------------------------------------------------------------------
    function get_option_value(self,key) result(val)
        class(bcproperty_t),    intent(in)  :: self
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
