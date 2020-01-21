!> Dictionary container
!!
!!  @author Nathan A. Wukie
!!  @date   2/1/2016
!!
!!   The general structure for this type was crafted after
!!   an example on fortranwiki.org, hash table example.
!!   However, here we are just traversing a linked-list to find things,
!!   and not hashing anything
!!
!------------------------------------------------------------------------
module type_dict
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    implicit none
    private






    !> Linked-List container type - used by Dictionary
    !!      - character:real pair
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !---------------------------------------------------
    type :: llreal_t
        ! Node pointer
        type(llreal_t), pointer :: child => null()

        ! Character-Real value pair
        character(len=:), allocatable   :: key
        real(rk)                        :: val


    contains
        procedure   :: set  => set_llreal
        procedure   :: get  => get_llreal
        procedure   :: free => free_llreal

    end type llreal_t
    !***************************************************





    !> Linked-List container type - used by Dictionary
    !!      - character:integer pair
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------
    type :: llint_t
        ! Node pointer
        type(llint_t), pointer :: child => null()

        ! Character-Integer value pair
        character(len=:), allocatable   :: key
        integer(ik)                     :: val


    contains
        procedure   :: set  => set_llint
        procedure   :: get  => get_llint
        procedure   :: free => free_llint

    end type llint_t
    !*****************************************************




    !> Linked-List container type - used by Dictionary
    !!      - character:logical pair
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   12/27/2018
    !!
    !------------------------------------------------------
    type :: lllogical_t
        ! Node pointer
        type(lllogical_t), pointer :: child => null()

        ! Character-Integer value pair
        character(len=:), allocatable   :: key
        logical                         :: val


    contains
        procedure   :: set  => set_lllogical
        procedure   :: get  => get_lllogical
        procedure   :: free => free_lllogical

    end type lllogical_t
    !*****************************************************




    !> Linked-List container type - used by Dictionary
    !!      - character:character pair
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   12/27/2018
    !!
    !------------------------------------------------------
    type :: llcharacter_t
        ! Node pointer
        type(llcharacter_t), pointer :: child => null()

        ! Character-Integer value pair
        character(:),   allocatable :: key
        character(:),   allocatable :: val


    contains
        procedure   :: set  => set_llcharacter
        procedure   :: get  => get_llcharacter
        procedure   :: free => free_llcharacter

    end type llcharacter_t
    !*****************************************************



    !> Dictionary Type for storing key-value pairs
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-------------------------------------------------
    type, public :: dict_t

        type(llreal_t)      :: llreal
        type(llint_t)       :: llint
        type(lllogical_t)   :: lllogical
        type(llcharacter_t) :: llcharacter

    contains
        procedure           :: print
        procedure           :: contains

        generic             :: set => set_real, set_int, set_logical, set_character
        generic             :: get => get_real, get_int, get_logical, get_character

        procedure, private  :: set_real => set_real_dict
        procedure, private  :: get_real => get_real_dict

        procedure, private  :: set_int  => set_int_dict
        procedure, private  :: get_int  => get_int_dict

        procedure, private  :: set_logical  => set_logical_dict
        procedure, private  :: get_logical  => get_logical_dict

        procedure, private  :: set_character  => set_character_dict
        procedure, private  :: get_character  => get_character_dict

    end type dict_t
    !*************************************************



contains










    !>  Print dictionary contents
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !-------------------------------------------------------------------------------
    subroutine print(self)
        class(dict_t),  intent(in), target  :: self

        type(llint_t),          pointer  :: llint  => null()
        type(llreal_t),         pointer  :: llreal => null()
        type(lllogical_t),      pointer  :: lllogical => null()
        type(llcharacter_t),    pointer  :: llcharacter => null()


    
        llint       => self%llint
        llreal      => self%llreal
        lllogical   => self%lllogical
        llcharacter => self%llcharacter

        ! If the current node contains an initialized key-value pair
        ! and it happens to be the key we are looking for
        do while ( associated(llint) ) 
            ! Print current node if allocated.
            if ( allocated(llint%key) ) then
                call write_line(llint%key,llint%val)
            end if

            ! Move to next node, if it exists.
            if ( associated(llint%child) ) then
                llint => llint%child
            else
                llint => null()
            end if
        end do


        ! If the current node contains an initialized key-value pair
        ! and it happens to be the key we are looking for
        do while ( associated(llreal) ) 
            ! Print current node if allocated.
            if ( allocated(llreal%key) ) then
                call write_line(llreal%key,llreal%val)
            end if

            ! Move to next node, if it exists.
            if ( associated(llreal%child) ) then
                llreal => llreal%child
            else
                llreal => null()
            end if
        end do


        ! If the current node contains an initialized key-value pair
        ! and it happens to be the key we are looking for
        do while ( associated(lllogical) ) 
            ! Print current node if allocated.
            if ( allocated(lllogical%key) ) then
                call write_line(lllogical%key,lllogical%val)
            end if

            ! Move to next node, if it exists.
            if ( associated(lllogical%child) ) then
                lllogical => lllogical%child
            else
                lllogical => null()
            end if
        end do


        ! If the current node contains an initialized key-value pair
        ! and it happens to be the key we are looking for
        do while ( associated(llcharacter) ) 
            ! Print current node if allocated.
            if ( allocated(llcharacter%key) ) then
                call write_line(llcharacter%key,llcharacter%val)
            end if

            ! Move to next node, if it exists.
            if ( associated(llcharacter%child) ) then
                llcharacter => llcharacter%child
            else
                llcharacter => null()
            end if
        end do


    end subroutine print
    !********************************************************************************







    !>  Return logical indicating if a key is registered in the dictionary.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !---------------------------------------------------------------------------------
    function contains(self,key) result(key_status)
        class(dict_t),  intent(in), target  :: self
        character(*),   intent(in)          :: key

        type(llint_t),      pointer  :: llint  => null()
        type(llreal_t),     pointer  :: llreal => null()
        type(lllogical_t),  pointer  :: lllogical => null()
        type(llcharacter_t), pointer  :: llcharacter => null()

        logical :: key_status, key_found


        llint       => self%llint
        llreal      => self%llreal
        lllogical   => self%lllogical
        llcharacter => self%llcharacter


        key_status = .false.    ! initialize to fail unless key is found

        ! Traverse nodes in linked-list
        do while ( associated(llint) )

            ! Check if key matches current node%key
            key_found  = .false.
            if ( allocated(llint%key) ) then
                key_found = ( key == llint%key )
            end if

            if ( key_found ) then
                key_status = .true.
                exit
            else
                llint => llint%child
            end if

        end do


        ! Traverse nodes in linked-list
        do while ( associated(llreal) )

            ! Check if key matches current node%key
            key_found  = .false.
            if ( allocated(llreal%key) ) then
                key_found = ( key == llreal%key )
            end if

            if ( key_found ) then
                key_status = .true.
                exit
            else
                llreal => llreal%child
            end if

        end do


        ! Traverse nodes in linked-list
        do while ( associated(lllogical) )

            ! Check if key matches current node%key
            key_found  = .false.
            if ( allocated(lllogical%key) ) then
                key_found = ( key == lllogical%key )
            end if

            if ( key_found ) then
                key_status = .true.
                exit
            else
                lllogical => lllogical%child
            end if

        end do


        ! Traverse nodes in linked-list
        do while ( associated(llcharacter) )

            ! Check if key matches current node%key
            key_found  = .false.
            if ( allocated(llcharacter%key) ) then
                key_found = ( key == llcharacter%key )
            end if

            if ( key_found ) then
                key_status = .true.
                exit
            else
                llcharacter => llcharacter%child
            end if

        end do


    end function contains
    !*********************************************************************************


















    !
    !   Routines for Dictionary
    !
    subroutine get_real_dict(self,key,val)
        class(dict_t),    intent(in)    :: self
        character(len=*), intent(in)    :: key
        real(kind=rk),    intent(inout) :: val

        call self%llreal%get(key,val)
    end subroutine 

    subroutine set_real_dict(self,key,val)
        class(dict_t),    intent(inout) :: self
        character(len=*), intent(in)    :: key
        real(kind=rk),    intent(in)    :: val

        call self%llreal%set(key,val)

    end subroutine




    subroutine get_int_dict(self,key,val)
        class(dict_t),     intent(in)    :: self
        character(len=*),  intent(in)    :: key
        integer(kind=ik),  intent(inout) :: val

        call self%llint%get(key,val)
    end subroutine 

    subroutine set_int_dict(self,key,val)
        class(dict_t),      intent(inout) :: self
        character(len=*),   intent(in)    :: key
        integer(kind=ik),   intent(in)    :: val

        call self%llint%set(key,val)

    end subroutine



    subroutine get_logical_dict(self,key,val)
        class(dict_t),    intent(in)    :: self
        character(len=*), intent(in)    :: key
        logical,          intent(inout) :: val

        call self%lllogical%get(key,val)
    end subroutine 

    subroutine set_logical_dict(self,key,val)
        class(dict_t),    intent(inout) :: self
        character(len=*), intent(in)    :: key
        logical,          intent(in)    :: val

        call self%lllogical%set(key,val)
    end subroutine





    subroutine get_character_dict(self,key,val)
        class(dict_t),      intent(in)    :: self
        character(*),       intent(in)    :: key
        character(*),       intent(inout) :: val

        call self%llcharacter%get(key,val)
    end subroutine get_character_dict

    subroutine set_character_dict(self,key,val)
        class(dict_t),  intent(inout) :: self
        character(*),   intent(in)    :: key
        character(*),   intent(in)    :: val

        call self%llcharacter%set(key,val)
    end subroutine set_character_dict




    !>  Set key/value pair for real
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/2/2016
    !!
    !---------------------------------------------------------------------------------------------
    recursive subroutine set_llreal(list,key,val)
        class(llreal_t),    intent(inout) :: list
        character(*),       intent(in)    :: key
        real(rk),           intent(in)    :: val

        ! Check if the current node is an allocated pair
        if ( allocated(list%key) ) then
            !check if the key matches what we are looking to set
            if ( list%key == key ) then
                list%val = val
            else
                ! The key does not match what we are looking for. Check for child node
                if ( .not. associated(list%child) ) then
                    allocate(list%child)    ! Create node for new key/val pair
                end if
                ! Now that we have allocated an empty node we can set it's properties
                call set_llreal(list%child,key,val)
            end if
        else
            ! So, the current node exists, but the key-value pair needs set
            list%key = key
            list%val = val
        end if

    end subroutine set_llreal
    !*****************************************************************************************



    !-----------------------------------------------------------------------------------------
    recursive subroutine get_llreal(list,key,val)
        class(llreal_t),    intent(in)      :: list
        character(*),       intent(in)      :: key
        real(kind=rk),      intent(inout)   :: val

        ! If the current node contains an initialized key-value pair
        ! and it happens to be the key we are looking for
        if (allocated(list%key) .and. (list%key == key)) then
            ! Return the associated value
            val = list%val
        ! The current node did not contain the key we were looking for.
        ! Check if there is another node to move to.
        else if (associated(list%child)) then
            !There was another node, so call get on that node
            call get_llreal(list%child,key,val)
        ! We searched the whole list and found no valid key.
        ! BAD!
        else
            call chidg_signal_one(FATAL,"dict%get: key was not found",key)
        end if

    end subroutine get_llreal
    !******************************************************************************************


    recursive subroutine free_llreal(list)
        class(llreal_t),    intent(inout)   :: list

        if (associated(list%child)) then
            call free_llreal(list%child)
            deallocate(list%child)
        end if
        list%child => null()
        if (allocated(list%key)) deallocate(list%key)

    end subroutine free_llreal



    !==============================================
    !
    !   character:logical   key:value
    !
    !==============================================



    !>  Set key/value pair for logical
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   12/27/2018
    !!
    !---------------------------------------------------------------------------------------------
    recursive subroutine set_lllogical(list,key,val)
        class(lllogical_t), intent(inout) :: list
        character(*),       intent(in)    :: key
        logical,            intent(in)    :: val

        ! Check if the current node is an allocated pair
        if ( allocated(list%key) ) then
            !check if the key matches what we are looking to set
            if ( list%key == key ) then
                list%val = val
            else
                ! The key does not match what we are looking for. Check for child node
                if ( .not. associated(list%child) ) then
                    allocate(list%child)    ! Create node for new key/val pair
                end if
                ! Now that we have allocated an empty node we can set it's properties
                call set_lllogical(list%child,key,val)
            end if
        else
            ! So, the current node exists, but the key-value pair needs set
            list%key = key
            list%val = val
        end if

    end subroutine set_lllogical
    !*****************************************************************************************



    !-----------------------------------------------------------------------------------------
    recursive subroutine get_lllogical(list,key,val)
        class(lllogical_t), intent(in)      :: list
        character(*),       intent(in)      :: key
        logical,            intent(inout)   :: val

        ! If the current node contains an initialized key-value pair
        ! and it happens to be the key we are looking for
        if (allocated(list%key) .and. (list%key == key)) then
            ! Return the associated value
            val = list%val
        ! The current node did not contain the key we were looking for.
        ! Check if there is another node to move to.
        else if (associated(list%child)) then
            !There was another node, so call get on that node
            call get_lllogical(list%child,key,val)
        else
            call chidg_signal_one(FATAL,"dict%get: key was not found",key)
        end if

    end subroutine get_lllogical
    !******************************************************************************************


    recursive subroutine free_lllogical(list)
        class(lllogical_t),    intent(inout)   :: list

        if (associated(list%child)) then
            call free_lllogical(list%child)
            deallocate(list%child)
        end if
        list%child => null()
        if (allocated(list%key)) deallocate(list%key)

    end subroutine free_lllogical



    !==============================================
    !
    !   character:character key:value
    !
    !==============================================



    !>  Set key/value pair for character 
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   12/27/2018
    !!
    !---------------------------------------------------------------------------------------------
    recursive subroutine set_llcharacter(list,key,val)
        class(llcharacter_t),   intent(inout) :: list
        character(*),           intent(in)    :: key
        character(*),           intent(in)    :: val

        ! Check if the current node is an allocated pair
        if ( allocated(list%key) ) then
            !check if the key matches what we are looking to set
            if ( list%key == key ) then
                list%val = val
            else
                ! The key does not match what we are looking for. Check for child node
                if ( .not. associated(list%child) ) then
                    allocate(list%child)    ! Create node for new key/val pair
                end if
                ! Now that we have allocated an empty node we can set it's properties
                call set_llcharacter(list%child,key,val)
            end if
        else
            ! So, the current node exists, but the key-value pair needs set
            list%key = key
            list%val = val
        end if

    end subroutine set_llcharacter
    !*****************************************************************************************




    !-----------------------------------------------------------------------------------------
    recursive subroutine get_llcharacter(list,key,val)
        class(llcharacter_t),   intent(in)      :: list
        character(*),           intent(in)      :: key
        character(*),           intent(inout)   :: val

        ! If the current node contains an initialized key-value pair
        ! and it happens to be the key we are looking for
        if (allocated(list%key) .and. (list%key == key)) then
            ! Return the associated value
            val = list%val

        ! The current node did not contain the key we were looking for.
        ! Check if there is another node to move to.
        else if (associated(list%child)) then
            !There was another node, so call get on that node
            call get_llcharacter(list%child,key,val)

        else
            call chidg_signal_one(FATAL,"dict%get: key was not found",key)
        end if

    end subroutine get_llcharacter
    !******************************************************************************************



    recursive subroutine free_llcharacter(list)
        class(llcharacter_t),    intent(inout)   :: list

        if (associated(list%child)) then
            call free_llcharacter(list%child)
            deallocate(list%child)
        end if
        list%child => null()
        if (allocated(list%key)) deallocate(list%key)

    end subroutine free_llcharacter












    !===============================================
    !
    !   Routines for character-integer linked list
    !
    !===============================================

    !> Find the existing key and set val, or create a new node
    !! and set the associated value.
    !!
    !!  @author Nathan A. Wukie
    !!
    !------------------------------------------------------------
    recursive subroutine set_llint(list,key,val)
        class(llint_t),      intent(inout) :: list
        character(len=*),    intent(in)    :: key
        integer(kind=ik),    intent(in)    :: val

        integer(kind=ik)                :: keylen


        keylen = len(key)

        ! Check if the current node is an allocated pair
        if ( allocated(list%key) ) then
            !check if the key matches what we are looking to set
            if ( list%key == key ) then
                list%val = val
            else
                ! The key does not match what we are looking for. Check for child node
                if ( .not. associated(list%child) ) then
                    allocate(list%child)    ! Create node for new key/val pair
                end if
                ! Now that we have allocated an empty node we can set it's properties
                call set_llint(list%child,key,val)
            end if
        else
            ! So, the current node exists, but the key-value pair needs set
            list%key = key
            list%val = val
        end if

    end subroutine set_llint
    !************************************************************
    


    !> Traverse the linked-list to find the 'key' and return the associated value
    !!
    !!
    !!  @author Nathan A. Wukie
    !!
    !----------------------------------------------------------------------------
    recursive subroutine get_llint(list,key,val)
        class(llint_t),     intent(in)      :: list
        character(len=*),   intent(in)      :: key
        integer(kind=ik),   intent(inout)   :: val

        ! If the current node contains an initialized key-value pair
        ! and it happens to be the key we are looking for
        if (allocated(list%key) .and. (list%key == key)) then
            ! Return the associated value
            val = list%val

        ! The current node did not contain the key we were looking for.
        ! Check if there is another node to move to.
        else if (associated(list%child)) then
            !There was another node, so call get on that node
            call get_llint(list%child,key,val)

        ! We searched the whole list and found no valid key.
        ! BAD!
        else
            call chidg_signal_one(FATAL,"dict%get: key was not found",key)
        end if

    end subroutine get_llint
    !****************************************************************************






    !> Remove a node from the list
    !!
    !!
    !----------------------------------------------------------------------------
    recursive subroutine free_llint(list)
        class(llint_t),    intent(inout)   :: list

        if (associated(list%child)) then
            call free_llint(list%child)
            deallocate(list%child)
        end if
        list%child => null()
        if (allocated(list%key)) deallocate(list%key)

    end subroutine


end module type_dict
