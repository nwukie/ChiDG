module type_datawrapper
    implicit none

    type, public :: datawrapper_t
        class(*), allocatable :: elem
    end type datawrapper_t

end module type_datawrapper



module type_vector
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use type_datawrapper,   only: datawrapper_t
    implicit none




    type, public :: vector_t
        integer(ik)                         :: size_        = 0
        integer(ik)                         :: capacity_    = 0
        integer(ik)                         :: buffer_      = 20




        type(datawrapper_t),   allocatable :: data(:)

    contains
        procedure, public   :: size
        procedure, public   :: capacity


        !< Data modifiers
        procedure,  public  :: push_back
        procedure,  public  :: clear
        procedure           :: increase_capacity


        !< Data accessors
        procedure, public   :: at
    end type vector_t



contains



    !> This function returns the number of elements stored in the container
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/17/2016
    !!
    !!
    !-------------------------------------------------------------------------
    function size(self) result(res)
        class(vector_t),    intent(in)  :: self

        integer(ik) :: res

        res = self%size_
    end function



    !> This function returns the total capacity of the container to store elements
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/17/2016
    !!
    !!
    !-------------------------------------------------------------------------
    function capacity(self) result(res)
        class(vector_t),    intent(in)  :: self

        integer(ik) :: res

        res = self%capacity_
    end function









    !> Store element to end of vector
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/17/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine push_back(self,element)
        class(vector_t), intent(inout)   :: self
        class(*),        intent(in)      :: element
        
        type(datawrapper_t)     :: wrapper
        logical                 :: capacity_reached
        integer                 :: ierr, size

        !
        ! Test if container has storage available. If not, then increase capacity
        !
        capacity_reached = (self%size() == self%capacity())
        if (capacity_reached) then
            call self%increase_capacity()
        end if
        

        !
        ! Allocate wrapper component and store data
        !
        !allocate(wrapper%elem, source=element, stat=ierr)
        !if (ierr /= 0) call AllocationError


        !
        ! Add element to end of vector
        !
        size = self%size()
        allocate(self%data(size + 1)%elem, source=element, stat=ierr)
        if (ierr /= 0) call AllocationError
        !self%data(size + 1) = element


        !
        ! Increment number of stored elements
        !
        self%size_ = self%size_ + 1


    end subroutine push_back
    !*****************************************************************************************









    !> Clear container contents
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !-------------------------------------------------------------------------
    subroutine clear(self)
        class(vector_t),   intent(inout)   :: self

        self%size_     = 0
        self%capacity_ = 0

        deallocate(self%data)

    end subroutine clear












    !> Access element at index location
    !!
    !!
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    function at(self,index) result(res)
        class(vector_t),    intent(in)  :: self
        integer(ik),        intent(in)  :: index

        class(*), allocatable   :: res
        logical                 :: out_of_bounds

        !
        ! Check vector bounds
        !
        out_of_bounds = (index > self%size())
        if (out_of_bounds) then
            call chidg_signal(FATAL,'vector_t%at: out of bounds access')
        end if


        !
        ! Allocate result
        !

        allocate(res, source=self%data(index)%elem)


    end function



















    !> Increase the storage capacity of the vector by a buffer size predefined in the container
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/17/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine increase_capacity(self)
        class(vector_t),    intent(inout)   :: self

        type(datawrapper_t), allocatable    :: temp(:)
        integer(ik)                         :: newsize, ierr


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
    !******************************************************************************************
















end module type_vector
