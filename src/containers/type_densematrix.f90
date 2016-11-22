!> Data type for storing the dense block matrices for the linearization of each element
!!  @author Nathan A. Wukie
module type_densematrix
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: NO_PROC
    implicit none




    !> Storage for full dense/full matrices.
    !!
    !!  block associativity
    !!  Domain-global index of the element, with which this block is associated.
    !!  For example, a given element has linearization blocks xi_min,xi_max,eta_min, etc.
    !!
    !!
    !!  Block dimensions
    !!   -->j
    !!  |
    !!  v
    !!  i
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/6/2016
    !!  @note   Extended for parallel
    !!
    !!
    !------------------------------------------------------------------------------------------
    type, public :: densematrix_t


        ! zero value indicates unassigned
        integer(ik)     :: parent_proc_ = NO_PROC
        integer(ik)     :: dparent_g_   = 0   !< Global domain index of the element matrix was linearized with respect to
        integer(ik)     :: dparent_l_   = 0   !< Local domain index of the element matrix was linearized with respect to
        integer(ik)     :: eparent_g_   = 0   !< Global element index of the element matrix was linearized with respect to
        integer(ik)     :: eparent_l_   = 0   !< Local element index of the element matrix was linearized with respect to

        ! imat index of transposed densematrix
        integer(ik)     :: itranspose_  = 0

        ! If associated parent data is being received from another processor, location in chidgVector%recv to find it
        integer(ik)     :: recv_comm    = 0
        integer(ik)     :: recv_domain  = 0
        integer(ik)     :: recv_element = 0

        ! Block storage
        integer(ik)             :: nrows_
        integer(ik)             :: ncols_
        real(rk),   allocatable :: mat(:,:)

    contains
        ! Initializers
        procedure :: init               !< Initialize block with general-sized matrix storage

        ! Getters
        procedure :: dparent_g          !< return parent domain
        procedure :: dparent_l          !< return parent domain
        procedure :: eparent_g          !< return parent element
        procedure :: eparent_l          !< return parent element
        procedure :: parent_proc        !< return the processor rank of the parent
        procedure :: itranspose         !< return imat index of transposed densematrix
        procedure :: nentries           !< return number of matrix entries
        procedure :: idim               !< return i-dimension of matrix storage
        procedure :: jdim               !< return j-dimension of matrix storage
        procedure :: dump               !< print out matrix contents
        procedure :: get_recv_comm
        procedure :: get_recv_domain
        procedure :: get_recv_element

        ! Setters
        procedure :: resize             !< resize matrix storage
        procedure :: set_recv_comm
        procedure :: set_recv_domain
        procedure :: set_recv_element
        procedure :: set_itranspose


    end type densematrix_t
    !******************************************************************************************





    private

contains





    !> Subroutine for initializing general dense-block storage
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!  @param[in]  idim        Column-dimension of dense matrix to be initialized
    !!  @param[in]  jbim        Row-dimension of dense matrix to be initialized
    !!  @param[in]  dparent     Integer index of parent domain
    !!  @param[in]  eparent     Integer index of parent element
    !!
    !------------------------------------------------------------------------------------------
    subroutine init(self,idim,jdim,dparent_g,dparent_l,eparent_g,eparent_l,parent_proc)
        class(densematrix_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: idim
        integer(ik),            intent(in)      :: jdim
        integer(ik),            intent(in)      :: dparent_g
        integer(ik),            intent(in)      :: dparent_l
        integer(ik),            intent(in)      :: eparent_g
        integer(ik),            intent(in)      :: eparent_l
        integer(ik),            intent(in)      :: parent_proc

        integer(ik) :: ierr

        !
        ! Set parents
        !
        self%dparent_g_   = dparent_g
        self%dparent_l_   = dparent_l
        self%eparent_g_   = eparent_g
        self%eparent_l_   = eparent_l
        self%parent_proc_ = parent_proc


        !
        ! Set matrix size
        !
        self%nrows_ = idim
        self%ncols_ = jdim

        !
        ! Allocate block storage
        ! Check if storage was already allocated and reallocate if necessary
        !
        if ( allocated(self%mat) ) deallocate(self%mat)
        allocate(self%mat(idim,jdim),stat=ierr)
        if (ierr /= 0) call AllocationError


        !
        ! Initialize to zero
        !
        self%mat = 0._rk

    end subroutine init
    !*******************************************************************************************








    !> return i-dimension of block storage
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!  
    !!  @return i   Integer of the column-dimension of the stored matrix
    !!
    !-------------------------------------------------------------------------------------------
    function idim(self) result(i)
        class(densematrix_t), intent(in)    :: self
        integer(ik)                         :: i

        i = size(self%mat,1)

    end function idim
    !*******************************************************************************************






    !> return j-dimension of block storage
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @return j   Integer of the row-dimension of the stored matrix
    !!
    !------------------------------------------------------------------------------------------
    function jdim(self) result(j)
        class(densematrix_t), intent(in)    :: self
        integer(ik)                         :: j

        j = size(self%mat,2)

    end function jdim
    !******************************************************************************************








    !> return number of entries in block storage
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @return     nentries    Integer of the total number of matrix entries stored
    !!
    !------------------------------------------------------------------------------------------
    function nentries(self) result(n)
        class(densematrix_t), intent(in)    :: self
        integer(ik)                         :: n

        n = size(self%mat,1) * size(self%mat,2)

    end function nentries
    !******************************************************************************************












    !> return index of parent domain
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @return     par     Integer index of the parent domain.
    !!
    !------------------------------------------------------------------------------------------
    function dparent_g(self) result(par)
        class(densematrix_t),   intent(in)  :: self
        integer(ik)                         :: par

        par = self%dparent_g_

    end function dparent_g
    !******************************************************************************************




    !> return index of parent domain
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @return     par     Integer index of the parent domain.
    !!
    !------------------------------------------------------------------------------------------
    function dparent_l(self) result(par)
        class(densematrix_t),   intent(in)  :: self
        integer(ik)                         :: par

        par = self%dparent_l_

    end function dparent_l
    !******************************************************************************************







    !> return index of parent element
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @return     par     Integer index of the parent element.
    !!
    !------------------------------------------------------------------------------------------
    function eparent_g(self) result(par)
        class(densematrix_t), intent(in) :: self
        integer(ik)                     :: par

        par = self%eparent_g_

    end function eparent_g
    !******************************************************************************************



    !> return index of parent element
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @return     par     Integer index of the parent element.
    !!
    !------------------------------------------------------------------------------------------
    function eparent_l(self) result(par)
        class(densematrix_t), intent(in) :: self
        integer(ik)                      :: par

        par = self%eparent_l_

    end function eparent_l
    !******************************************************************************************





    !> return index of parent element
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @return     par     Integer index of the parent element.
    !!
    !------------------------------------------------------------------------------------------
    function parent_proc(self) result(par)
        class(densematrix_t), intent(in) :: self
        integer(ik)                      :: par

        par = self%parent_proc_

    end function parent_proc
    !******************************************************************************************






    !>  Return itranspose, imat index of transposed densematrix.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/22/2016
    !!
    !!  @return     itranspose     Integer index(imat) of the densematrix in the transposed location.
    !!
    !------------------------------------------------------------------------------------------
    function itranspose(self) result(par)
        class(densematrix_t), intent(in) :: self
        integer(ik)                      :: par

        par = self%itranspose_

    end function itranspose
    !******************************************************************************************









    !> Resize dense-block storage
    !!
    !! @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine resize(self,idim,jdim)
        class(densematrix_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: idim
        integer(ik),            intent(in)      :: jdim

        integer(ik) :: ierr


        !
        ! Set matrix size
        !
        self%nrows_ = idim
        self%ncols_ = jdim


        !
        ! Allocate block storage
        ! Check if storage was already allocated and reallocate if necessary
        !
        if (allocated(self%mat)) then
            deallocate(self%mat)
            allocate(self%mat(idim,jdim),stat=ierr)
        else
            allocate(self%mat(idim,jdim),stat=ierr)
        end if
        if (ierr /= 0) call AllocationError

    end subroutine resize
    !******************************************************************************************







    !>  Set recv_comm component.
    !!
    !!  recv_comm is the comm index in a chidgVector that information coming from the parent
    !!  element can be found.
    !!
    !!  chidgVector%recv%comm(recv_comm)%dom(recv_domain)%vec(recv_element)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine set_recv_comm(self,recv_comm)
        class(densematrix_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: recv_comm

        self%recv_comm = recv_comm

    end subroutine set_recv_comm
    !******************************************************************************************




    !>  Set recv_domain component.
    !!
    !!  recv_domain is the comm index in a chidgVector that information coming from the parent
    !!  element can be found.
    !!
    !!  chidgVector%recv%comm(recv_comm)%dom(recv_domain)%vec(recv_element)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine set_recv_domain(self,recv_domain)
        class(densematrix_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: recv_domain

        self%recv_domain = recv_domain

    end subroutine set_recv_domain
    !******************************************************************************************





    !>  Set recv_element component.
    !!
    !!  recv_element is the comm index in a chidgVector that information coming from the parent
    !!  element can be found.
    !!
    !!  chidgVector%recv%comm(recv_comm)%dom(recv_domain)%vec(recv_element)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine set_recv_element(self,recv_element)
        class(densematrix_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: recv_element

        self%recv_element = recv_element

    end subroutine set_recv_element
    !******************************************************************************************





    !>  Get recv_comm component.
    !!
    !!  recv_comm is the comm index in a chidgVector that information coming from the parent
    !!  element can be found.
    !!
    !!  chidgVector%recv%comm(recv_comm)%dom(recv_domain)%vec(recv_element)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_recv_comm(self) result(recv_comm)
        class(densematrix_t),   intent(inout)   :: self

        integer(ik) :: recv_comm

        recv_comm = self%recv_comm

    end function get_recv_comm
    !******************************************************************************************






    !>  Get recv_domain component.
    !!
    !!  recv_domain is the domain index in a chidgVector that information coming from the parent
    !!  element can be found.
    !!
    !!  chidgVector%recv%comm(recv_comm)%dom(recv_domain)%vec(recv_element)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_recv_domain(self) result(recv_domain)
        class(densematrix_t),   intent(inout)   :: self

        integer(ik) :: recv_domain

        recv_domain = self%recv_domain

    end function get_recv_domain
    !******************************************************************************************






    !>  Get recv_element component.
    !!
    !!  recv_element is the element index in a chidgVector that information coming from the parent
    !!  element can be found.
    !!
    !!  chidgVector%recv%comm(recv_comm)%dom(recv_domain)%vec(recv_element)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/15/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function get_recv_element(self) result(recv_element)
        class(densematrix_t),   intent(inout)   :: self

        integer(ik) :: recv_element

        recv_element = self%recv_element

    end function get_recv_element
    !******************************************************************************************






    !>  Set itranspose component.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/22/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine set_itranspose(self,itranspose)
        class(densematrix_t),   intent(inout)   :: self
        integer(ik),            intent(in)      :: itranspose

        self%itranspose_ = itranspose

    end subroutine set_itranspose
    !******************************************************************************************







    !>  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !-------------------------------------------------
    subroutine dump(self)
        class(densematrix_t), intent(inout) :: self


        integer(ik) :: irow

        print*, self%dparent_g_
        print*, self%dparent_l_
        print*, self%eparent_g_
        print*, self%eparent_l_

        do irow = 1,size(self%mat,1)
            print*, self%mat(irow,:)
        end do


    end subroutine dump
    !*************************************************












end module type_densematrix
