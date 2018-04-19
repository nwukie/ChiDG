module type_densematrix_vector
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: ZERO
    use type_densematrix,   only: densematrix_t
    use type_element,       only: element_t
    use mod_chidg_mpi,      only: IRANK
    use DNAD_D
    implicit none


    !>  Vector container for storing dynamically sized arrays of integers.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    type, public :: densematrix_vector_t

        integer(ik)                 :: idomain_g
        integer(ik)                 :: idomain_l
        integer(ik)                 :: ielement_g
        integer(ik)                 :: ielement_l

        integer(ik)                 :: size_      = 0
        integer(ik)                 :: capacity_  = 0
        integer(ik)                 :: buffer_    = 7

        type(densematrix_t),    allocatable :: data_(:)
        real (rk),              allocatable :: mass(:,:)

    contains
        ! Initializers
        procedure, public   :: init                 ! initialize with size and capacity = 0


        procedure, public   :: capacity             ! return the current allocated capacity
        procedure, public   :: loc                  ! return the location of a stored value
        procedure, public   :: size => data_size    ! return the number of stored elements.
                                                    ! Associate to avoid clash with intrinsic.


        ! Data modifiers
        procedure, public   :: push_back
        procedure, public   :: clear
        procedure, public   :: setzero
        procedure, public   :: store
        procedure, public   :: store_element
        procedure, private  :: increase_capacity

        procedure, public   :: set_itranspose
        procedure, public   :: set_recv_comm
        procedure, public   :: set_recv_domain
        procedure, public   :: set_recv_element
        procedure, public   :: get_mass_matrix
        procedure, public   :: get_recv_comm
        procedure, public   :: get_recv_domain
        procedure, public   :: get_recv_element


        ! Data accessors
        procedure, public   :: at               ! return data from element densematrix_vector%at(ielem)
        procedure, public   :: data             ! return full data vector
        procedure, public   :: dmat             ! return densematrix array from element densematrix_vector%dmat(ielem)
        procedure, public   :: dparent_g        ! return parent global domain for the index position densematrix
        procedure, public   :: eparent_g        ! return parent global element for the index position densematrix
        procedure, public   :: dparent_l        ! return parent local domain for the index position densematrix
        procedure, public   :: eparent_l        ! return parent local element for the index position densematrix
        procedure, public   :: tparent          ! return parent time level for the index position densematrix
        procedure, public   :: parent_proc      ! return parent processor rank
        procedure, public   :: itranspose       ! return itranspose, imat index of densematrix in transposed location.
        procedure, public   :: get_diagonal     ! return index of densematrix representing the diagonal.
        procedure, public   :: get_idomain_g    ! return the global domain index which the densematrix_vector is associated with
        procedure, public   :: get_ielement_g   ! return the global element index which the densematrix vector is associated with
        procedure, public   :: get_idomain_l    ! return the local domain index which the densematrix_vector is associated with
        procedure, public   :: get_ielement_l   ! return the local element indexwhich the densematrix vector is associated with


        ! Process
        procedure, public   :: restrict

    end type densematrix_vector_t
    !*****************************************************************************************



contains


    !>  Subroutine for initializing general densematrix_vector
    !!
    !!  @author Matteo Ugolotti
    !!  @date 11/26/2016
    !!
    !!  @param[in]  idomain_g       global domain index which the densematrix_vector is associated with
    !!  @param[in]  ielement_g      global element index which the densematrix_vector is associated with
    !!  @param[in]  idomain_l       local domain index which the densematrix_vector is associated with
    !!  @param[in]  ielement_l      local element index which the densematrix_vector is associated with
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine init(self,idomain_g,ielement_g,idomain_l,ielement_l)
        class(densematrix_vector_t),    intent(inout)   :: self
        integer(ik),                    intent(in)      :: idomain_g,ielement_g,idomain_l,ielement_l
        
        !
        ! set domain and element index
        !

        self%idomain_g  = idomain_g
        self%ielement_g = ielement_g
        self%idomain_l  = idomain_l
        self%ielement_l = ielement_l

        !
        ! set size and capacity to zero
        !
        
        self%size_      = ZERO
        self%capacity_  = ZERO


    end subroutine init
    !*****************************************************************************************





    !>  This function returns idomain_g for the densematrix_vector
    !!
    !!  @author Matteo Ugolotti
    !!  @date 11/26/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_idomain_g(self) result(res)
        class(densematrix_vector_t),    intent(inout)   :: self
        integer(ik)      :: res
        
        res = self%idomain_g

    end function get_idomain_g
    !*****************************************************************************************






    !>  This function returns idomain_l for the densematrix_vector
    !!
    !!  @author Matteo Ugolotti
    !!  @date 11/26/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_idomain_l(self) result(res)
        class(densematrix_vector_t),    intent(inout)   :: self
        integer(ik)      :: res
        
        res = self%idomain_l

    end function get_idomain_l
    !*****************************************************************************************






    !>  This function returns ielement_g for the densematrix_vector
    !!
    !!  @author Matteo Ugolotti
    !!  @date 11/26/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_ielement_g(self) result(res)
        class(densematrix_vector_t),    intent(inout)   :: self
        integer(ik)      :: res
        
        res = self%ielement_g

    end function get_ielement_g
    !*****************************************************************************************






    !>  This function returns ielement_l for the densematrix_vector
    !!
    !!  @author Matteo Ugolotti
    !!  @date 11/26/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_ielement_l(self) result(res)
        class(densematrix_vector_t),    intent(inout)   :: self
        integer(ik)      :: res
        
        res = self%ielement_l

    end function get_ielement_l
    !*****************************************************************************************






    !>  This function returns the number of elements stored in the container. 
    !!
    !!  Normally, the name here would just be 'size', but we were getting conflicts with 
    !!  the intrinsic function 'size', so it was renamed data_size. Then, it is bound
    !!  to the procedure size in the data type. size => data_size.
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti + Nathan A. Wukie
    !!  @date   11/07/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function data_size(self) result(res)
        class(densematrix_vector_t),   intent(in)  :: self

        integer(ik) :: res

        res = self%size_

    end function data_size
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
    !------------------------------------------------------------------------------------------
    function loc(self,idomain_g,ielem_g,itime) result(res)
        class(densematrix_vector_t),    intent(in)  :: self
        integer(ik),                    intent(in)  :: idomain_g
        integer(ik),                    intent(in)  :: ielem_g
        integer(ik),                    intent(in)  :: itime

        type(densematrix_t) :: temp_mat
        integer(ik)         :: res,ival

        res = 0
        do ival = 1,self%size()
            temp_mat = self%at(ival)
            if ( temp_mat%dparent_g() == idomain_g .and. &
                 temp_mat%eparent_g() == ielem_g   .and. &
                 temp_mat%tparent() == itime ) then 
                res = ival
                exit
            end if
        end do

    end function loc
    !******************************************************************************************





    !> Store element at end of vector
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti + Nathan A. Wukie
    !!  @date   11/07/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine push_back(self,element)
        class(densematrix_vector_t),   intent(inout)   :: self
        type(densematrix_t),           intent(in)      :: element

        logical     :: capacity_reached
        integer(ik) :: size
        integer(ik) :: duplicate

        ! Test whether the densematrix we are pushing is already there
        duplicate = self%loc(element%dparent_g(),element%eparent_g(),element%tparent())

        if ( duplicate == 0 ) then
              
            ! Test if container has storage available. If not, then increase capacity
            capacity_reached = (self%size() == self%capacity())
            if (capacity_reached) then
                call self%increase_capacity()
            end if

            ! Add element to end of vector
            size = self%size()
            self%data_(size + 1) = element

            ! Increment number of stored elements
            self%size_ = self%size_ + 1

        else 
           call chidg_signal_two(FATAL,"densematrix_vector_t%push_back: attempt to push the same densematrix twice,operation not allowed.", element%dparent_g(), element%eparent_g())
        end if


    end subroutine push_back
    !*******************************************************************************************







    !> Clear container contents
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti + Nathan A. Wukie
    !!  @date   11/07/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine clear(self)
        class(densematrix_vector_t),   intent(inout)   :: self

        self%size_     = 0
        self%capacity_ = 0
        
        if (allocated(self%data_)) deallocate(self%data_)

    end subroutine clear
    !*******************************************************************************************






    !> Set all densematrix_vector_t storage to zero
    !!
    !!  @author Matteo Ugolotti
    !!  @date 11/14/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
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
    !******************************************************************************************






    !> Access element at index location
    !!
    !!  @author Mayank Sharma + Matteo Ugolotti + Nathan A. Wukie
    !!  @date   11/07/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    function at(self,index) result(res)
        class(densematrix_vector_t),    intent(in)  :: self
        integer(ik),                    intent(in)  :: index

        logical             :: out_of_bounds
        type(densematrix_t) :: res

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
        
        real(rk),   allocatable     :: res(:,:)
        logical                     :: out_of_bounds
        
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



    !>  Set the mass matrix attribute of the densematrix vector equal to the mass matrix
    !!  coming in from an object of type_element
    !!
    !!  @Mayank Sharma + Matteo Ugolotti
    !!  @date   2/13/2017
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine get_mass_matrix(self,elem)
        class(densematrix_vector_t),    intent(inout)   :: self
        class(element_t),               intent(in)      :: elem


        self%mass = elem%mass


    end subroutine get_mass_matrix
    !****************************************************************************************




    !>  Return dparent_g() function from densematrix
    !!
    !!  @Matteo Ugolotti
    !!  @date   11/12/2016
    !!
    !----------------------------------------------------------------------------------------
    function dparent_g(self,index) result (par)
        class(densematrix_vector_t),    intent(in)  :: self
        integer(ik),                    intent(in)  :: index

        integer(ik) :: par

        par = self%data_(index)%dparent_g()

    end function dparent_g
    !****************************************************************************************







    !>  Return eparent_g() function from a densematrix
    !!
    !!  @Matteo Ugolotti
    !!  @date   11/12/2016
    !!
    !----------------------------------------------------------------------------------------
    function eparent_g(self,index) result (par)
        class(densematrix_vector_t),    intent(in)  :: self
        integer(ik),                    intent(in)  :: index

        integer(ik) :: par
        
        par = self%data_(index)%eparent_g()

    end function eparent_g
    !****************************************************************************************










    !>  Return dparent_l() function from densematrix
    !!
    !!  @Matteo Ugolotti
    !!  @date   11/12/2016
    !!
    !----------------------------------------------------------------------------------------
    function dparent_l(self,index) result (par)
        class(densematrix_vector_t),    intent(in)  :: self
        integer(ik),                    intent(in)  :: index

        integer(ik) :: par

        par = self%data_(index)%dparent_l()

    end function dparent_l
    !****************************************************************************************







    !>  Return eparent_l() function from a densematrix
    !!
    !!  @Matteo Ugolotti
    !!  @date   11/12/2016
    !!
    !----------------------------------------------------------------------------------------
    function eparent_l(self,index) result (par)
        class(densematrix_vector_t),    intent(in)  :: self
        integer(ik),                    intent(in)  :: index

        integer(ik) :: par
        
        par = self%data_(index)%eparent_l()

    end function eparent_l
    !****************************************************************************************



    !>  Return tparent() function from a densematrix
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/18/2018
    !!
    !----------------------------------------------------------------------------------------
    function tparent(self,index) result (par)
        class(densematrix_vector_t),    intent(in)  :: self
        integer(ik),                    intent(in)  :: index

        integer(ik) :: par
        
        par = self%data_(index)%tparent()

    end function tparent
    !****************************************************************************************






    !>  Return parent processor index for a densematrix.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !----------------------------------------------------------------------------------------
    function parent_proc(self,index) result(proc)
        class(densematrix_vector_t),    intent(in)  :: self
        integer(ik),                    intent(in)  :: index

        integer(ik) :: proc

        proc = self%data_(index)%parent_proc()

    end function parent_proc
    !****************************************************************************************





    !>  Return itranspose for a densematrix.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !----------------------------------------------------------------------------------------
    function itranspose(self,index) result(proc)
        class(densematrix_vector_t),    intent(in)  :: self
        integer(ik),                    intent(in)  :: index

        integer(ik) :: proc

        proc = self%data_(index)%itranspose()

    end function itranspose
    !****************************************************************************************




    !>  Set itranspose for a given densematrix.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine set_itranspose(self,index,itranspose)
        class(densematrix_vector_t),    intent(inout)   :: self
        integer(ik),                    intent(in)      :: index
        integer(ik),                    intent(in)      :: itranspose

        call self%data_(index)%set_itranspose(itranspose)

    end subroutine set_itranspose
    !****************************************************************************************






    !>  Set recv_comm component for a given densematrix.
    !!
    !!  recv_comm is the comm index in a chidgVector that information coming from the parent
    !!  element can be found.
    !!
    !!  chidgVector%recv%comm(recv_comm)%dom(recv_domain)%vec(recv_element)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine set_recv_comm(self,index,recv_comm)
        class(densematrix_vector_t),    intent(inout)   :: self
        integer(ik),                    intent(in)      :: index
        integer(ik),                    intent(in)      :: recv_comm

        call self%data_(index)%set_recv_comm(recv_comm)

    end subroutine set_recv_comm
    !****************************************************************************************



    !>  Set recv_domain component for a given densematrix.
    !!
    !!  recv_domain is the domain index in a chidgVector that information coming from the parent
    !!  element can be found.
    !!
    !!  chidgVector%recv%comm(recv_comm)%dom(recv_domain)%vec(recv_element)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine set_recv_domain(self,index,recv_domain)
        class(densematrix_vector_t),    intent(inout)   :: self
        integer(ik),                    intent(in)      :: index
        integer(ik),                    intent(in)      :: recv_domain

        call self%data_(index)%set_recv_domain(recv_domain)

    end subroutine set_recv_domain
    !****************************************************************************************




    !>  Set recv_element component for a given densematrix.
    !!
    !!  recv_element is the domain index in a chidgVector that information coming from the parent
    !!  element can be found.
    !!
    !!  chidgVector%recv%comm(recv_comm)%dom(recv_domain)%vec(recv_element)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !---------------------------------------------------------------------------------------
    subroutine set_recv_element(self,index,recv_element)
        class(densematrix_vector_t),    intent(inout)   :: self
        integer(ik),                    intent(in)      :: index
        integer(ik),                    intent(in)      :: recv_element

        call self%data_(index)%set_recv_element(recv_element)

    end subroutine set_recv_element
    !****************************************************************************************




    !>  Get recv_comm component for a given densematrix.
    !!
    !!  recv_comm is the comm index in a chidgVector that information coming from the parent
    !!  element can be found.
    !!
    !!  chidgVector%recv%comm(recv_comm)%dom(recv_domain)%vec(recv_element)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !---------------------------------------------------------------------------------------
    function get_recv_comm(self,index) result(recv_comm)
        class(densematrix_vector_t),    intent(inout)   :: self
        integer(ik),                    intent(in)      :: index

        integer(ik) :: recv_comm

        recv_comm = self%data_(index)%get_recv_comm()

    end function get_recv_comm
    !****************************************************************************************






    !>  Get recv_domain component for a given densematrix.
    !!
    !!  recv_domain is the domain index in a chidgVector that information coming from the parent
    !!  element can be found.
    !!
    !!  chidgVector%recv%comm(recv_comm)%dom(recv_domain)%vec(recv_element)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !---------------------------------------------------------------------------------------
    function get_recv_domain(self,index) result(recv_domain)
        class(densematrix_vector_t),    intent(inout)   :: self
        integer(ik),                    intent(in)      :: index

        integer(ik) :: recv_domain

        recv_domain = self%data_(index)%get_recv_domain()

    end function get_recv_domain
    !****************************************************************************************







    !>  Get recv_element component for a given densematrix.
    !!
    !!  recv_element is the element index in a chidgVector that information coming from the parent
    !!  element can be found.
    !!
    !!  chidgVector%recv%comm(recv_comm)%dom(recv_domain)%vec(recv_element)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !---------------------------------------------------------------------------------------
    function get_recv_element(self,index) result(recv_element)
        class(densematrix_vector_t),    intent(inout)   :: self
        integer(ik),                    intent(in)      :: index

        integer(ik) :: recv_element

        recv_element = self%data_(index)%get_recv_element()

    end function get_recv_element
    !****************************************************************************************









    !>  Return the index of the densematrix corresponding to the diagonal.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_diagonal(self) result(index)
        class(densematrix_vector_t),    intent(in)   :: self

        integer(ik) :: imat, index


        index = 0
        do imat = 1,self%size()
            if ( (self%data_(imat)%dparent_g() == self%idomain_g) .and. &
                 (self%data_(imat)%eparent_g() == self%ielement_g) ) then
                 index = imat
                 exit
            end if
        end do

    end function get_diagonal
    !*****************************************************************************************









    !> Store derivative data to the densematrix at index location in the densematrix_vector_t
    !!
    !!  @param[in]  vector  Array of modes from the spatial scheme, with embedded partial 
    !!                      derivatives for the linearization matrix
    !!
    !!  @Matteo Ugolotti
    !!  @date   11/12/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine store(self,index,ivar,nterms,integral)
        class(densematrix_vector_t),    intent(inout)   :: self
        integer(ik),                    intent(in)      :: index
        integer(ik),                    intent(in)      :: ivar
        integer(ik),                    intent(in)      :: nterms
        type(AD_D),                     intent(in)      :: integral(:)

        integer(ik) :: iarray,irow,irow_start


        ! Compute correct row offest for ivar
        irow_start = ( (ivar-1) * nterms )
        
        ! Loop through integral values, for each value store its derivative.
        ! The integral values here should be components of the RHS vector. 
        ! An array of partial derivatives from an AD_A variable should be stored 
        ! as a row in the block matrix
        do iarray = 1,size(integral)
            ! Do a += operation to add derivaties to any that are currently stored
            irow = irow_start + iarray
            self%data_(index)%mat(irow,:) = self%data_(index)%mat(irow,:) + integral(iarray)%xp_ad_
        end do

    end subroutine store
    !****************************************************************************************




    !>  Direct matrix contribution instead of pulling derivatives out of AD
    !!  data structures.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/19/2018
    !!
    !------------------------------------------------------------------------------
    subroutine store_element(self,index,contribution)
        class(densematrix_vector_t),    intent(inout)   :: self
        integer(ik),                    intent(in)      :: index
        real(rk),                       intent(in)      :: contribution(:,:)

        self%data_(index)%mat = self%data_(index)%mat + contribution

    end subroutine store_element
    !*******************************************************************************





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
        
        ! Get number of stored elements
        size = self%size()

        ! Allocate result
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







    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/23/2017
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function restrict(self,nterms_r) result(restricted)
        class(densematrix_vector_t),    intent(in)  :: self
        integer(ik),                    intent(in)  :: nterms_r

        type(densematrix_vector_t)  :: restricted
        type(densematrix_t)         :: dmatrix
        integer(ik)                 :: imat


        !
        ! Restrict mass matrix
        !
        restricted%mass = self%mass(1:nterms_r,1:nterms_r)


        !
        ! Store element identifier
        !
        restricted%idomain_g  = self%idomain_g
        restricted%idomain_l  = self%idomain_l
        restricted%ielement_g = self%ielement_g
        restricted%ielement_l = self%ielement_l


        !
        ! Loop through self%data_, for each block, restrict and add to new restricted object
        !
        do imat = 1,self%size()
            ! Get matrix
            dmatrix = self%at(imat)

            ! Push back restricted densematrix onto restricted densematrix_vector
            call restricted%push_back(dmatrix%restrict(nterms_r))
        end do



    end function restrict
    !*****************************************************************************************









end module type_densematrix_vector
