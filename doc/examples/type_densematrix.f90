!> Data type for storing the dense block matrices for the linearization of each element
!!  @author Nathan A. Wukie
module type_densematrix
#include <messenger.h>
    use mod_kinds,  only: rk,ik
    implicit none




    !> Storage for full dense/full matrices.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !--------------------------------------------------------------------------------------------------------
    !> [densematrix_t]
    type, public :: densematrix_t
        ! ZERO VALUE INDICATES UNASSIGNED
        integer(ik), private    :: dparent_ = 0     !< parent domain    For chimera reference accross domains
        integer(ik), private    :: eparent_ = 0     !< parent element

        ! Block storage
        ! NOTE: Assumes square blocks. TODO: extend for variable element solution expansion.
        real(rk),  dimension(:,:), allocatable :: mat


    end type densematrix_t
    !> [densematrix_t]
    !*************************************************************************************************************





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
    !------------------------------------------------------------------------------------------------------------------
    subroutine init_gen(self,idim,jdim,dparent,eparent)
        class(densematrix_t), intent(inout)  :: self
        integer(ik),         intent(in)     :: idim, jdim, dparent, eparent

        integer(ik) :: ierr

        !
        ! Set parents
        !
        self%dparent_ = dparent
        self%eparent_ = eparent



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


        !
        ! Initialize to zero
        !
        self%mat = 0._rk

    end subroutine
    !******************************************************************************************************************










    !> Subroutine for initializing square dense-block storage
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!  @param[in]  bsize       Row and column dimension for square matrix to be initialized
    !!  @param[in]  dparent     Integer index of parent domain
    !!  @param[in]  eparent     Integer index of parent element
    !!
    !--------------------------------------------------------------------------------------------------------------------
    subroutine init_square(self,bsize,dparent,eparent)
        class(densematrix_t), intent(inout)  :: self
        integer(ik),         intent(in)     :: bsize, dparent, eparent

        integer(ik) :: ierr



        !
        ! Block parents
        !
        self%dparent_ = dparent
        self%eparent_ = eparent



        !
        ! Allocate block storage
        ! Check if storage was already allocated and reallocate if necessary
        !
        if (allocated(self%mat)) then
            deallocate(self%mat)
            allocate(self%mat(bsize,bsize),stat=ierr)
        else
            allocate(self%mat(bsize,bsize),stat=ierr)
        end if
        if (ierr /= 0) call AllocationError



        !
        ! Initialize to zero
        !
        self%mat = 0._rk

    end subroutine
    !*******************************************************************************************************************












    !> return i-dimension of block storage
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!  
    !!  @return i   Integer of the column-dimension of the stored matrix
    !!
    !-------------------------------------------------------------------------------------------------------------------
    function idim(self) result(i)
        class(densematrix_t), intent(in)    :: self
        integer(ik)                         :: i

        i = size(self%mat,1)

    end function
    !*******************************************************************************************************************













    !> return j-dimension of block storage
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @return j   Integer of the row-dimension of the stored matrix
    !!
    !-------------------------------------------------------------------------------------------------------------------
    function jdim(self) result(j)
        class(densematrix_t), intent(in)    :: self
        integer(ik)                         :: j

        j = size(self%mat,2)

    end function
    !*******************************************************************************************************************













    !> return number of entries in block storage
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @return     nentries    Integer of the total number of matrix entries stored
    !!
    !-------------------------------------------------------------------------------------------------------------------
    function nentries(self) result(n)
        class(densematrix_t), intent(in)    :: self
        integer(ik)                         :: n

        n = size(self%mat,1) * size(self%mat,2)

    end function
    !*******************************************************************************************************************












    !> return index of parent domain
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @return     par     Integer index of the parent domain.
    !!
    !-------------------------------------------------------------------------------------------------------------------
    function dparent(self) result(par)
        class(densematrix_t),   intent(in)  :: self
        integer(ik)                         :: par

        par = self%dparent_

    end function
    !*******************************************************************************************************************











    !> return index of parent element
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @return     par     Integer index of the parent element.
    !!
    !-------------------------------------------------------------------------------------------------------------------
    function eparent(self) result(par)
        class(densematrix_t), intent(in) :: self
        integer(ik)                     :: par

        par = self%eparent_
    end function
    !*******************************************************************************************************************














    !> Resize dense-block storage
    !!
    !! @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-------------------------------------------------------------------------------------------------------------------
    subroutine resize(self,idim,jdim)
        class(densematrix_t), intent(inout)  :: self
        integer(ik),         intent(in)     :: idim,jdim

        integer(ik) :: ierr

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

    end subroutine
    !*******************************************************************************************************************







    !> set index of parent
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !-------------------------------------------------------------------------------------------------------------------
    subroutine reparent(self,par)
        class(densematrix_t), intent(inout)  :: self
        integer(ik),         intent(in)     :: par

        self%eparent_ = par

    end subroutine
    !*******************************************************************************************************************










    !>  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !-------------------------------------------------
    subroutine dump(self)
        class(densematrix_t), intent(inout) :: self


        integer(ik) :: irow

        print*, self%dparent_
        print*, self%eparent_

        do irow = 1,size(self%mat,1)
            print*, self%mat(irow,:)
        end do


    end subroutine dump
    !*************************************************












end module type_densematrix
