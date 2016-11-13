module type_densematrix_vector
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use type_densematrix,   only: densematrix_t
    use type_densevector,   only: densevector_t
    implicit none

    !TODO: Add helper routine dparent() to get index
    !TODO: Add helper routine store to help with storage in blockmatrix
    !TODO: Add helper routine clear to help clear in blockmatrix


    !>  Vector container for storing dynamically sized arrays of integers.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    type, public :: densematrix_vector_t

        integer(ik)                 :: size_        = 0
        integer(ik)                 :: capacity_    = 0
        integer(ik)                 :: buffer_      = 7
        type(densematrix_t),   allocatable  :: data_(:)

    contains
        procedure, public   :: size     !< return the number of stored elements
        procedure, public   :: capacity !< return the current allocated capacity
        procedure, public   :: loc      !< return the location of a stored value


        ! Data modifiers
        procedure, public   :: push_back
        procedure, public   :: clear
        procedure, private  :: increase_capacity


        ! Data accessors
        procedure, public   :: at       !< return data from element densematrix_vector%at(ielem)
        procedure, public   :: data     !< return full data vector

    end type densematrix_vector_t
    !*****************************************************************************************



contains



    !> This function returns the number of elements stored in the container
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti + Nathan A. Wukie
    !!  @date   11/07/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function size(self) result(res)
        class(densematrix_vector_t),   intent(in)  :: self

        integer(ik) :: res

        res = self%size_
    end function size
    !*****************************************************************************************






    !> This function returns the total capacity of the container to store elements
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti + Nathan A. Wukie
    !!  @date 11/07/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function capacity(self) result(res)
        class(densematrix_vector_t),   intent(in)  :: self

        integer(ik) :: res

        res = self%capacity_
    end function capacity
    !*****************************************************************************************
















    !>  This function returns the location of a given value. If not found, returns 0
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti + Nathan A. Wukie
    !!  @date   11/07/2016
    !!
    !-------------------------------------------------------------------------------------------
    function loc(self,idomain_g,ielem_g) result(res)
        class(densematrix_vector_t),   intent(in)  :: self
        integer(ik),                   intent(in)  :: idomain_g,ielem_g

        type(densematrix_t)                        :: temp_mat
        integer(ik) :: res,ival

        res = 0


        do ival = 1,self%size()

            temp_mat = self%at(ival)
            
            if ( temp_mat%dparent_g() == idomain_g .and. temp_mat%eparent_g() == ielem_g ) then 
                res = ival
                exit
            end if

        end do

    end function loc
    !*******************************************************************************************





    !> Store element at end of vector
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti + Nathan A. Wukie
    !!  @date   11/07/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine push_back(self,element)
        class(densematrix_vector_t),   intent(inout)   :: self
        type(densematrix_t),           intent(in)      :: element

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














    !> Clear container contents
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti + Nathan A. Wukie
    !!  @date   11/07/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(densematrix_vector_t),   intent(inout)   :: self

        self%size_     = 0
        self%capacity_ = 0

        
        if (allocated(self%data_)) deallocate(self%data_)

    end subroutine clear
    !*********************************************************************************************











    !> Access element at index location
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti + Nathan A. Wukie
    !!  @date   11/07/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function at(self,index) result(res)
        class(densematrix_vector_t),   intent(in)  :: self
        integer(ik),        intent(in)  :: index

        type(densematrix_t) :: res
        logical             :: out_of_bounds

        !
        ! Check vector bounds
        !
        out_of_bounds = (index > self%size())
        if (out_of_bounds) then
            call chidg_signal(FATAL,"densematrix_vector_t%at: out of bounds access")
        end if


        !
        ! Allocate result
        !
        res = self%data_(index)

    end function at
    !****************************************************************************************





    !> Access entire data vector
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti + Nathan A. Wukie
    !!  @date   11/07/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function data(self) result(res)
        class(densematrix_vector_t),   intent(in)  :: self
        
        type(densematrix_t), allocatable    :: res(:)
        integer(ik)                         :: size
        
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
    !!  @author Mayank Sharma + Matteo Ugolotti + Nathan A. Wukie
    !!  @date   11/07/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine increase_capacity(self)
        class(densematrix_vector_t),   intent(inout)   :: self

        type(densematrix_t), allocatable    :: temp(:)
        integer(ik)                         :: newsize, ierr


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
















end module type_densematrix_vector
