module type_densematrix_vector
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mode_constants,     only: ZERO
    use type_densematrix,   only: densematrix_t
    use type_densevector,   only: densevector_t
    use DNAD
    implicit none


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
        procedure, public   :: setzero
        procedure, private  :: increase_capacity
        procedure, public   :: store_dmv


        ! Data accessors
        procedure, public   :: at               !< return data from element densematrix_vector%at(ielem)
        procedure, public   :: data             !< return full data vector
        procedure, public   :: dmat             !< return densematrix array from element densematrix_vector%dmat(ielem)
        procedure, public   :: dparent_g_dmvi   !< return parent domain for the index position densematrix
        procedure, public   :: eparent_g_dmv    !< return parent element for the index position densematrix
        procedure, public   :: find             !< find element in densematrix vector based idonor_domain_g ad i_element_g

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






    !> Set all densematrix_vector_t storage to zero
    !!
    !!  @author Matteo Ugolotti
    !!  @date 11/14/2016
    !!
    !!
    !--------------------------------------------------------------------------------------------
    subroutine setzero(self)
        class(densematrix_vector_t),    intent(inout)   :: self
         
        integer(ik)       :: ival
        

        do ival = 1, self%size()
            
            !
            ! Check if the densematrix is actually allocated
            !
                    
            if (allocated(self%data_(ival)%mat)) then

            self%data_(ival)%mat = ZERO

            end if

        end do


    end subroutine setzero
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







    !> Access densematrix at index location
    !!
    !!  @Matteo Ugolotti
    !!  @date   11/12/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function dmat(self,index) result (res)
        class(densematrix_vector_t),    intent(in)  :: self
        integer(ik),                    intent(in)  :: index
        
        real(rk),allocatable            :: res(:,:)
        logical                         :: out_of_bounds
        
        !
        ! Check vector bounds
        !
        out_of_bounds = (index > self%size())
        if (out_of_bounds) then
            call chidg_signal(FATAL,"densematrix_vector_t%dmat: out of bounds access")
        end if


        res = self%data_(index)%mat


    end function dmat
    !****************************************************************************************








    !> Inherit dparent_g() function from densematrix
    !!
    !!  @Matteo Ugolotti
    !!  @date   11/12/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function dparent_g_dmv(self,index) result (par)
        class(densematrix_vector_t),    intent(in)  :: self
        integer(ik),                    intent(in)  :: index
        integer(ik)                                 :: par
        


        par = self%data_(index)%dparent_g()


    end function dparent_g_dmv
    !****************************************************************************************







    !> Inherit eparent_g() function from densematrix
    !!
    !!  @Matteo Ugolotti
    !!  @date   11/12/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function eparent_g_dmv(self,index) result (par)
        class(densematrix_vector_t),    intent(in)  :: self
        integer(ik),                    intent(in)  :: index
        integer(ik)                                 :: par
        


        par = self%data_(index)%eparent_g()


    end function eparent_g_dmv
    !****************************************************************************************







    !> Find element in densematrix_vector based idonor_domain_g ad i_element_g
    !!
    !!  @Matteo Ugolotti
    !!  @date   11/14/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    function find(self,donor_domain,donor_element) result (res)
        class(densematrix_vector_t),    intent(in)  :: self
        integer(ik),                    intent(in)  :: donor_domain, donor_element
        real(rk),                       intent(out) :: res
        
        integer(ik)         :: ival
        logical             :: matrix_match = .false.
        logical             :: no_donor_matrix = .false. 

        res = 0

        do ival = 1 , self%size()

            matrix_match = ( (donor_domain == self%dparent_g_dmv(ival)) .and. &
                            (donor_element == self%eparent_g_dmv(ival)) )

            if ( matrix_match ) then
                res = ival
                exit
            end if
            

        end do ! ival
        
    
        no_donor_matrix = (res == 0)
        if (no_donor_matrix) call chidg_signal(MSG, 'densematrix_vector%find: no donor densematrix found to store derivative')


        

    end function find
    !****************************************************************************************






    !> Store derivative data to the densematrix at index location in the densematrix_vector_t
    !!
    !!  @param[in]  vector  Array of modes from the spatial scheme, withembedded partial derivatives for the linearization matrix
    !!
    !!
    !!
    !!
    !!
    !!  @Matteo Ugolotti
    !!  @date   11/12/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine store_dmv(self,index,ivar,nterms,vector)
        class(densematrix_vector_t),    intent(inout)   :: self
        type(AD_D),                     intent(in)      :: vector(:)
        integer(ik),                    intent(in)      :: ivar,nterms,index
        
        integer(ik)             :: iarray,irow_start


        !
        ! Compute correct row offest for ivar
        !

        irow_start = ( (ivar -1) * nterms)
        
        !
        ! Loop through integral values, for each value store its derivative.
        ! The integral values here should be components of the RHS vector. An array of partial derivatives from an AD_A variable
        ! should be stored as a row in the block matrix


        do iarray = 1,size(vector)

            !
            ! Do a += operation to add derivaties to any that are currently stored

            irow = irow_start + iarray

            self%data_(index)%mat(irow,:) = self%data_(index)%mat(irow,:) + vector(iarray)%xp_ad_
        
        end do


    end subroutine store_dmv
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
