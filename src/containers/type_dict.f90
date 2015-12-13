!> Dictionary container
!!
!!   The general structure for this type was crafted after
!!   an example on fortranwiki.org, hash table example.
!!   However, here we are just traversing a linked-list to find things,
!!   and not hashing anything
!!
!------------------------------------------------------------------------
module type_dict
    use mod_kinds,      only: rk,ik
    implicit none
    private


    !> Linked-List container type - used by Dictionary
    !!      - character:real pair
    !!
    !!  @author Nathan A. Wukie
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





    !> Linked-List container type - used by Dictionary
    !!      - character:integer pair
    !!
    !!  @author Nathan A. Wukie
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





    !> Dictionary Type for storing key-value pairs
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !-------------------------------------------------
    type, public :: dict_t

        type(llreal_t)   :: llreal
        type(llint_t)   :: llint

    contains
        generic :: set => set_real, set_int
        generic :: get => get_real, get_int

        procedure, private   :: set_real => set_real_dict
        procedure, private   :: get_real => get_real_dict

        procedure, private   :: set_int  => set_int_dict
        procedure, private   :: get_int  => get_int_dict

    end type dict_t

!    interface set
!        module procedure set_real, set_int
!    end interface


contains

    !
    !   Routines for Dictionary
    !
    subroutine get_real_dict(self,key,val)
        class(dict_t),    intent(inout) :: self
        character(len=*), intent(in)    :: key
        real(kind=rk),    intent(out)   :: val

        call self%llreal%get(key,val)
    end subroutine 

    subroutine set_real_dict(self,key,val)
        class(dict_t),    intent(inout) :: self
        character(len=*), intent(in)    :: key
        real(kind=rk),    intent(in)    :: val

        call self%llreal%set(key,val)

    end subroutine




    subroutine get_int_dict(self,key,val)
        class(dict_t),     intent(inout) :: self
        character(len=*),  intent(in)    :: key
        integer(kind=ik),  intent(out)   :: val

        call self%llint%get(key,val)
    end subroutine 

    subroutine set_int_dict(self,key,val)
        class(dict_t),      intent(inout) :: self
        character(len=*),   intent(in)    :: key
        integer(kind=ik),   intent(in)    :: val

        call self%llint%set(key,val)

    end subroutine





    !===============================================
    !
    !   Routines for character-real linked list
    !
    !===============================================
    recursive subroutine set_llreal(list,key,val)
        class(llreal_t),     intent(inout) :: list
        character(len=*),    intent(in)    :: key
        real(kind=rk),       intent(in)    :: val

        integer(kind=ik)                :: keylen


        keylen = len(key)

        ! Check if the current node is an allocated pair
        if (allocated(list%key)) then
            ! If it is, check if the key matches what we are looking to set
            if (list%key /= key) then
                ! The key does not match what we are looking for,
                ! so check if there is another node
                if (.not. associated(list%child)) then
                    ! There was not another node, so we need to create one
                    ! in order to set the value
                    allocate(list%child)
                end if

                ! Now that we have allocated an empty node we can set it's properties
                call set_llreal(list%child,key,val)

            end if
        else
            if (.not. allocated(list%key)) then
                allocate (character(len=keylen) :: list%key)
            end if
        ! So, the current node exists, but the key-value pair hasn't bet allocated
!            allocate(character(len=keylen) :: list%key)
            list%key = key
            list%val = val
        end if

    end subroutine



    recursive subroutine get_llreal(list,key,val)
        class(llreal_t),    intent(in)  :: list
        character(len=*),   intent(in)  :: key
        real(kind=rk),      intent(out) :: val

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
            print*, "Error dict_t: key ", key, "was not found."
            stop
        end if

    end subroutine


    recursive subroutine free_llreal(list)
        class(llreal_t),    intent(inout)   :: list

        if (associated(list%child)) then
            call free_llreal(list%child)
            deallocate(list%child)
        end if
        list%child => null()
        if (allocated(list%key)) deallocate(list%key)

    end subroutine







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
        if (allocated(list%key)) then
            ! If it is, check if the key matches what we are looking to set
            if (list%key /= key) then
                ! The key does not match what we are looking for,
                ! so check if there is another node
                if (.not. associated(list%child)) then
                    ! There was not another node, so we need to create one
                    ! in order to set the value
                    allocate(list%child)
                end if

                ! Now that we have allocated an empty node we can set it's properties
                call set_llint(list%child,key,val)

            end if
        else
        ! So, the current node exists, but the key-value pair hasn't bet allocated
            if (.not. allocated(list%key)) then
                allocate (character(len=keylen) :: list%key)
            end if
            list%key = key
            list%val = val
        end if

    end subroutine
    


    !> Traverse the linked-list to find the 'key' and return the associated value
    !!
    !!
    !!  @author Nathan A. Wukie
    !!
    !----------------------------------------------------------------------------
    recursive subroutine get_llint(list,key,val)
        class(llint_t),     intent(in)  :: list
        character(len=*),   intent(in)  :: key
        integer(kind=ik),   intent(out) :: val

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
            print*, "Error dict_t: key ", key, "was not found."
            stop
        end if

    end subroutine






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
