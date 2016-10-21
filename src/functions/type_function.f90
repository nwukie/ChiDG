module type_function
#include <messenger.h>
    use mod_kinds,  only: rk, ik
    use type_dict,  only: dict_t
    implicit none
    private


    !>  Function class.
    !!
    !!  f = f(t,\vec{x}) 
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------
    type, public, abstract :: function_t

        character(len=:),       allocatable, private    :: name_
        character(len=1024),    allocatable             :: option_keys(:)
        real(rk),               allocatable             :: option_vals(:)

    contains

        ! Function initialization
        procedure(init_interface),      deferred    :: init             !< Initialize function and register options
        procedure                                   :: add_option       !< Add option to the function definition. Only use for initialization

        procedure                                   :: set_name         !< Add function name. Only use for initialization
        procedure                                   :: get_name         !< Return the function name

        ! Function set
        procedure                                   :: set_option       !< Set option value

        ! Function get
        procedure                                   :: get_noptions     !< Return number of options in the function
        procedure                                   :: get_option_key   !< Return the name of an option,  given an index
        procedure                                   :: get_option_value !< Return the value of an option, given a key

        ! Function check
        procedure                                   :: check_key_exists !< Return logical indicating if a key exists in the list.
        procedure                                   :: get_key_index    !< Return the index of a key in the current list.

        ! Function copute
        procedure(compute_interface),   deferred    :: compute          !< Elemental function definition

    end type function_t
    !********************************************************************************************





    abstract interface
        impure elemental function compute_interface(self,time,coord)
            use type_point, only: point_t
            use mod_kinds,  only: rk
            import function_t

            class(function_t),  intent(inout)   :: self
            real(rk),           intent(in)      :: time
            type(point_t),      intent(in)      :: coord
            real(rk)                            :: compute_interface
        end function
    end interface




    abstract interface
        subroutine init_interface(self)
            import function_t

            class(function_t),  intent(inout)  :: self

        end subroutine
    end interface




contains



    !>  Add a name for the function
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine set_name(self,fname)
        class(function_t),  intent(inout)   :: self
        character(*),       intent(in)      :: fname

        self%name_ = fname

    end subroutine set_name
    !********************************************************************************************







    !>  Return the name of the function.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !-----------------------------------------------------------------------------------------
    function get_name(self) result(fname)
        class(function_t),  intent(in)  :: self

        character(len=:),   allocatable :: fname

        fname = self%name_

    end function get_name
    !*****************************************************************************************







    !>  Procedure for adding function options to the container list. These can be set then 
    !!  via the 'set_option' procedure.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !!  @param[in]  key     String indicating a dictionary entry to set
    !!  @parma[in]  val     Real value to be associated with the key
    !!
    !--------------------------------------------------------------------------------------------
    subroutine add_option(self,key,val)
        class(function_t),  intent(inout)    :: self
        character(*),       intent(in)       :: key
        real(rk),           intent(in)       :: val
        
        logical     :: key_exists
        real(rk)    :: test
        integer(ik) :: noptions, iopt, key_location, ierr

        character(len=1024),    allocatable :: temp_keys(:)
        real(rk),               allocatable :: temp_vals(:)


        noptions = self%get_noptions()

        !
        ! Check if key already exists
        !
        key_exists = .false.
        if ( noptions > 0 ) then
            do iopt = 1,noptions
                
                !
                ! Test if current option key matches the incoming key
                !
                if ( trim(self%option_keys(iopt)) == trim(key) ) then
                    key_exists = .true.
                    key_location = iopt
                end if

            end do
        end if



        !
        ! Handle key_exists status
        !
        if ( key_exists ) then

            !
            ! Set value of existing key to incoming new value
            !
            self%option_vals(key_location) = val

        else

            !
            ! Extend the option allocation by 1 and set a new key/val pair.
            !
            allocate(temp_keys(noptions+1), stat=ierr)
            if (ierr /= 0) call AllocationError

            allocate(temp_vals(noptions+1), stat=ierr)
            if (ierr /= 0) call AllocationError


            !
            ! Copy existing key/val pairs
            !
            temp_keys(1:noptions) = self%option_keys(1:noptions)
            temp_vals(1:noptions) = self%option_vals(1:noptions)


            !
            ! Set new key/val pair
            !
            temp_keys(noptions+1) = key
            temp_vals(noptions+1) = val


            !
            ! Move allocation to self
            !
            call move_alloc(temp_keys, self%option_keys)
            call move_alloc(temp_vals, self%option_vals)

        end if



    end subroutine add_option
    !**********************************************************************************************









    !>  Set an option that exists in the available list of options. Option must have been previously
    !!  added with the 'add_option' procedure.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------
    subroutine set_option(self,key,val)
        class(function_t),  intent(inout)   :: self
        character(*),       intent(in)      :: key
        real(rk),           intent(in)      :: val

        integer(ik) :: iopt, noptions, key_loc
        logical     :: key_exists

        !
        ! Check that key exists in the current list.
        !
        key_exists = self%check_key_exists(key)


        !
        ! Set dictionary key/value pair
        !
        if (key_exists) then
            
            key_loc = self%get_key_index(key)
            self%option_vals(key_loc) = val

        else


            call write_line("Invalid function option '"//key//"' in "//self%get_name()//"%set('"//key//"',",val,')', delimiter='')
            call write_line(' ')
            call write_line('Valid function options are:')
            call write_line(' ')
            call write_line('key','current value')
            call write_line('----------------------------------------------')
            
            noptions = self%get_noptions()
            do iopt = 1,noptions
                call write_line(self%option_keys(iopt), self%option_vals(iopt))
            end do ! iopt

            call chidg_signal_one(FATAL,"function%set: key is not a valid option for the current function",key)

        end if

    end subroutine set_option
    !*********************************************************************************************
















    !>  Return an option key, given the index of an available option.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !---------------------------------------------------------------------------------------------
    function get_option_key(self,iopt) result(key)
        class(function_t),  intent(in)  :: self
        integer(ik),        intent(in)  :: iopt

        integer(ik)                     :: noptions
        character(len=:),   allocatable :: key

        noptions = self%get_noptions()

        if (iopt > noptions) then
            call chidg_signal(FATAL,"function%get_optionkey: option index exceeds the bounds of the available number of options.")
        else
            key = trim(self%option_keys(iopt))
        end if


    end function get_option_key
    !*********************************************************************************************












    !>  Return the value of an option, given an option key.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------
    function get_option_value(self,key) result(val)
        class(function_t),  intent(in)  :: self
        character(*),       intent(in)  :: key

        real(rk)        :: val
        integer(ik)     :: key_loc, noptions
        logical         :: key_exists

        noptions = self%get_noptions()

        key_exists = self%check_key_exists(key)



        !
        ! Handle key status
        !
        if ( key_exists ) then

            key_loc = self%get_key_index(key)

            val = self%option_vals(key_loc)

        else
            call chidg_signal(FATAL,"function%get_option_value: key not found")
        end if


    end function get_option_value
    !********************************************************************************************










    !>  Return the number of available options.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function get_noptions(self) result(nopt)
        class(function_t),  intent(in)  :: self

        integer(ik)     :: nopt


        if ( allocated(self%option_vals) ) then

            nopt = size(self%option_vals)
        else
            nopt = 0
        end if


    end function get_noptions
    !*********************************************************************************************













    !>  Given a key, check if it exists as a valid option for the function.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function check_key_exists(self,key) result(key_exists)
        class(function_t),  intent(in)  :: self
        character(*),       intent(in)  :: key

        logical     :: key_exists
        integer(ik) :: noptions, iopt


        noptions = self%get_noptions()


        !
        ! Check if key exists in the current list
        !
        key_exists = .false.
        if ( noptions > 0 ) then
            do iopt = 1,noptions
                
                !
                ! Test if current option key matches the incoming key
                !
                if ( trim(self%option_keys(iopt)) == trim(key) ) then
                    key_exists = .true.
                end if

            end do
        end if



    end function check_key_exists
    !*********************************************************************************************








    !>  Return the index of an option in the option list, given an identifying key.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/4/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------
    function get_key_index(self,key) result(key_loc)
        class(function_t),      intent(in)  :: self
        character(*),           intent(in)  :: key

        integer(ik)     :: noptions, iopt, key_loc

        noptions = self%get_noptions()


        key_loc = 0
        do iopt = 1,noptions

            if ( trim(key) == trim(self%option_keys(iopt)) ) then
                key_loc = iopt
                exit
            end if

        end do


        if ( key_loc == 0 ) then
            call chidg_signal_one(FATAL,"function%get_key_index: key not found in list", key)
        end if


    end function get_key_index
    !*********************************************************************************************











end module type_function
