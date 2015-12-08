module type_ivector
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    implicit none




    type, public :: ivector_t
        integer(ik)             :: size_        = 0
        integer(ik)             :: capacity_    = 0
        integer(ik)             :: buffer_      = 20




        integer(ik),   allocatable :: data_(:)

    contains
        procedure, public   :: size     ! return the number of stored elements
        procedure, public   :: capacity ! return the current allocated capacity
        procedure, public   :: loc      ! return the location of a stored value


        !< Data modifiers
        procedure, public   :: push_back
        procedure, public   :: clear
        procedure, private  :: increase_capacity


        !< Data accessors
        procedure, public   :: at       !< return data from element ivector%at(ielem)
        procedure, public   :: data     !< return full data vector

    end type ivector_t



contains



    !> This function returns the number of elements stored in the container
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !-------------------------------------------------------------------------
    function size(self) result(res)
        class(ivector_t),   intent(in)  :: self

        integer(ik) :: res

        res = self%size_
    end function



    !> This function returns the total capacity of the container to store elements
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !-------------------------------------------------------------------------
    function capacity(self) result(res)
        class(ivector_t),   intent(in)  :: self

        integer(ik) :: res

        res = self%capacity_
    end function
















    !> This function returns the location of a given value. If not found, returns 0
    !!
    !!  @author Nathan A. Wukie
    !!
    !--------------------------------------------------------------------------
    function loc(self,val) result(res)
        class(ivector_t),   intent(in)  :: self
        integer(ik),        intent(in)  :: val

        integer(ik) :: res, ival

        res = 0
        do ival = 1,self%size()

            if ( self%at(ival) == val ) then
                res = ival
                exit
            end if

        end do

    end function





    !> Store element at end of vector
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine push_back(self,element)
        class(ivector_t),   intent(inout)   :: self
        integer(ik),        intent(in)      :: element

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









    !> Clear container contents
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !-------------------------------------------------------------------------
    subroutine clear(self)
        class(ivector_t),   intent(inout)   :: self

        self%size_     = 0
        self%capacity_ = 0

        deallocate(self%data_)

    end subroutine clear










    !> Access element at index location
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !----------------------------------------------------------------------------------------
    function at(self,index) result(res)
        class(ivector_t),   intent(in)  :: self
        integer(ik),        intent(in)  :: index

        integer(ik) :: res
        logical     :: out_of_bounds

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

    end function





    !> Access entire data vector
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function data(self) result(res)
        class(ivector_t),   intent(in)  :: self
        
        integer(ik), allocatable    :: res(:)
        integer(ik)                 :: size
        
        !
        ! Get number of stored elements
        !
        size = self%size()

        !
        ! Allocate result
        !
        res = self%data_(1:size)


    end function













    !> Increase the storage capacity of the vector by a buffer size predefined in the container
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine increase_capacity(self)
        class(ivector_t),   intent(inout)   :: self

        integer(ik), allocatable    :: temp(:)
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
















end module type_ivector
