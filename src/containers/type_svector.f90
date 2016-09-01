module type_svector
#include <messenger.h>
    use mod_kinds,  only: rk, ik
    use mod_string, only: string_t
    implicit none




    !>  Vector container for storing dynamically sized arrays of integers.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    type, public :: svector_t

        integer(ik)                 :: size_        = 0
        integer(ik)                 :: capacity_    = 0
        integer(ik)                 :: buffer_      = 20
        type(string_t), allocatable :: data_(:)

    contains
        procedure, public   :: size     !< return the number of stored elements
        procedure, public   :: capacity !< return the current allocated capacity
        procedure, public   :: loc      !< return the location of a stored value


        ! Data modifiers
        procedure, public   :: push_back
        procedure, public   :: push_back_unique
        procedure, public   :: clear
        procedure, private  :: increase_capacity


        ! Data accessors
        procedure, public   :: at       !< return data from element svector%at(ielem)
        procedure, public   :: data     !< return full data vector

    end type svector_t
    !*****************************************************************************************



contains



    !> This function returns the number of elements stored in the container
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function size(self) result(res)
        class(svector_t),   intent(in)  :: self

        integer(ik) :: res

        res = self%size_
    end function size
    !*****************************************************************************************






    !> This function returns the total capacity of the container to store elements
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function capacity(self) result(res)
        class(svector_t),   intent(in)  :: self

        integer(ik) :: res

        res = self%capacity_
    end function capacity
    !*****************************************************************************************
















    !>  This function returns the location of a given value. If not found, returns 0
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-------------------------------------------------------------------------------------------
    function loc(self,string) result(res)
        class(svector_t),   intent(in)  :: self
        character(len=*),   intent(in)  :: string

        integer(ik)     :: res, istr
        type(string_t)  :: temp_string

        res = 0
        do istr = 1,self%size()

            ! Get current string
            temp_string = self%at(istr)

            if ( temp_string%str == string ) then
                res = istr
                exit
            end if

        end do

    end function loc
    !*******************************************************************************************





    !> Store element at end of vector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine push_back(self,element)
        class(svector_t),   intent(inout)   :: self
        type(string_t),     intent(in)      :: element

        logical     :: capacity_reached
        integer(ik) :: size


        !
        ! Test if container has storage available. If not, then increase capacity
        !
        capacity_reached = (self%size() == self%capacity())
        if (capacity_reached) then
            call self%increase_capacity()
        end if


        !
        ! Add element to end of vector
        !
        size = self%size()
        self%data_(size + 1) = element


        !
        ! Increment number of stored elements
        !
        self%size_ = self%size_ + 1


    end subroutine push_back
    !********************************************************************************************








    !>  Store element at end of vector, if it hasn't already been stored. No duplicate values.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   7/22/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine push_back_unique(self,element)
        class(svector_t),   intent(inout)   :: self
        type(string_t),     intent(in)      :: element

        integer(ik) :: loc
        logical     :: already_added


        !
        ! Check if element was already added
        !
        loc = self%loc(element%str)
        already_added = (loc /= 0)


        !
        ! If not already in the list, push back
        !
        if ( .not. already_added ) call self%push_back(element)

    end subroutine push_back_unique
    !********************************************************************************************





    !> Clear container contents
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(svector_t),   intent(inout)   :: self

        self%size_     = 0
        self%capacity_ = 0

        
        if (allocated(self%data_)) deallocate(self%data_)

    end subroutine clear
    !*********************************************************************************************











    !> Access element at index location
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function at(self,index) result(res)
        class(svector_t),   intent(in)  :: self
        integer(ik),        intent(in)  :: index

        type(string_t)  :: res
        logical         :: out_of_bounds

        !
        ! Check vector bounds
        !
        out_of_bounds = (index > self%size())
        if (out_of_bounds) then
            call chidg_signal(FATAL,"vector_t%at: out of bounds access")
        end if


        !
        ! Allocate result
        !
        res = self%data_(index)

    end function at
    !****************************************************************************************





    !> Access entire data vector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function data(self) result(res)
        class(svector_t),   intent(in)  :: self
        
        type(string_t), allocatable :: res(:)
        integer(ik)                 :: size
        
        !
        ! Get number of stored elements
        !
        size = self%size()

        !
        ! Allocate result
        !
        res = self%data_(1:size)


    end function data
    !*****************************************************************************************













    !> Increase the storage capacity of the vector by a buffer size predefined in the container
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine increase_capacity(self)
        class(svector_t),   intent(inout)   :: self

        type(string_t), allocatable :: temp(:)
        integer(ik)                 :: newsize, ierr


        !
        ! Allocate temporary vector of current size plus a buffer
        !
        if ( allocated(self%data_) ) then
            newsize = ubound(self%data_,1) + self%buffer_
        else
            newsize = self%buffer_
        end if

        allocate(temp(newsize),stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Copy any current data to temporary vector
        !
        if (allocated(self%data_)) then
            temp(lbound(self%data_,1):ubound(self%data_,1))  =  self%data_
        end if


        !
        ! Move alloc to move data back to self%data and deallocate temp
        !
        call move_alloc(FROM=temp,TO=self%data_)


        !
        ! Reset capacity info
        !
        self%capacity_ = newsize


    end subroutine increase_capacity
    !*****************************************************************************************







end module type_svector
