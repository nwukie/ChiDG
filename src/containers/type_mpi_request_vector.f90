module type_mpi_request_vector
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use mod_string,                 only: string_to_upper

    use mpi_f08,                    only: MPI_REQUEST
    implicit none



    !>  A vector class for storing a dynamic array of bc_t instances.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    type, public :: mpi_request_vector_t
        integer(ik)             :: size_        = 0
        integer(ik)             :: capacity_    = 0
        integer(ik)             :: buffer_      = 20

        type(mpi_request), allocatable :: data(:)


    contains

        procedure, public   :: size
        procedure, public   :: capacity

        ! Data modifiers
        procedure, public   :: push_back
        procedure, public   :: clear
        procedure, private  :: increase_capacity

        ! Data accessors
        procedure, public   :: at                   !< Return an instance from the specified index.

    end type mpi_request_vector_t
    !***************************************************************************************



contains



    !> This function returns the number of elements stored in the container
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    function size(self) result(res)
        class(mpi_request_vector_t),   intent(in)  :: self

        integer(ik) :: res

        res = self%size_
    end function size
    !**************************************************************************************







    !> This function returns the total capacity of the container to store elements
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    function capacity(self) result(res)
        class(mpi_request_vector_t),   intent(in)  :: self

        integer(ik) :: res

        res = self%capacity_
    end function capacity
    !**************************************************************************************








    !> Store element at end of vector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !--------------------------------------------------------------------------------------
    subroutine push_back(self,element)
        class(mpi_request_vector_t),    intent(inout)   :: self
        type(mpi_request),              intent(in)      :: element

        logical     :: capacity_reached
        integer(ik) :: size, ierr


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
        self%data(size + 1) = element


        !
        ! Increment number of stored elements
        !
        self%size_ = self%size_ + 1


    end subroutine push_back
    !***************************************************************************************








    !> Clear container contents
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !---------------------------------------------------------------------------------------
    subroutine clear(self)
        class(mpi_request_vector_t),   intent(inout)   :: self

        self%size_      = 0
        self%capacity_  = 0

        if (allocated(self%data)) deallocate(self%data)

    end subroutine clear
    !****************************************************************************************









    !> Access element at index location
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function at(self,index) result(res)
        class(mpi_request_vector_t),    intent(in)  :: self
        integer(ik),                    intent(in)  :: index

        integer             :: ierr
        type(mpi_request)   :: res
        logical             :: out_of_bounds

        !
        ! Check vector bounds
        !
        out_of_bounds = (index > self%size())
        if (out_of_bounds) then
            call chidg_signal(FATAL,'mpi_request_vector_t%at: out of bounds access')
        end if


        !
        ! Allocate result
        !
        res = self%data(index)

    end function at
    !*****************************************************************************************













    !> Increase the storage capacity of the vector by a buffer size predefined in the container
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine increase_capacity(self)
        class(mpi_request_vector_t),   intent(inout)   :: self

        type(mpi_request), allocatable   :: temp(:)
        integer(ik)                      :: newsize, ierr


        !
        ! Allocate temporary vector of current size plus a buffer
        !
        if ( allocated(self%data) ) then
            newsize = ubound(self%data,1) + self%buffer_
        else
            newsize = self%buffer_
        end if

        allocate(temp(newsize),stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Copy any current data to temporary vector
        !
        if (allocated(self%data)) then
            temp(lbound(self%data,1):ubound(self%data,1))  =  self%data
        end if


        !
        ! Move alloc to move data back to self%data and deallocate temp
        !
        call move_alloc(FROM=temp,TO=self%data)


        !
        ! Reset capacity info
        !
        self%capacity_ = newsize


    end subroutine increase_capacity
    !*******************************************************************************************







end module type_mpi_request_vector
