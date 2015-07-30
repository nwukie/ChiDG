!> Data type for storing the dense block matrices for the linearization of each element
!!  @author Nathan A. Wukie
module type_denseblock
    use mod_kinds,  only: rk,ik
    use messenger
    implicit none


    type, public :: denseblock_t
        !> block associativity
        !! Domain-global index of the element, with which this block is associated.
        !! For example, a given element has linearization blocks xi_min,xi_max,eta_min, etc.
        !! Below, we consider the linearization of element #6.
        !!
        !!   elem #5         elem #6         elem #7
        !!
        !!  blk_ximin       blk_diag        blk_ximax
        !!
        !!  The value of parent_ for the blocks would be:
        !! - blk_ximin = 5
        !! - blk_diag  = 6
        !! - blk_ximax = 7
        !!
        !! ZERO VALUE INDICATES UNASSIGNED
        integer(ik), private    :: parent_ = 0 !> parent element

        !> Block storage
        !! NOTE: Assumes square blocks, since this type is specifically
        real(rk),  dimension(:,:), allocatable :: mat

    contains
        !> Initializers
        generic, public :: init => init_gen, init_square
        procedure, private :: init_gen      !> Initialize block with general-sized matrix storage
        procedure, private :: init_square   !> Initialize block with square-sized matrix storage

        !> Getters
        !! Block dimensions
        !!
        !!      ---> j
        !!  |
        !!  |
        !!  v
        !!
        !!  i
        procedure :: parent     !> return parent element
        procedure :: nentries   !> return number of matrix entries
        procedure :: idim       !> return i-dimension of matrix storage
        procedure :: jdim       !> return j-dimension of matrix storage

        !> Setters
        procedure :: resize     !> resize matrix storage
        procedure :: reparent   !> reassign parent

        final :: destructor
    end type denseblock_t



    private
contains

    !> Subroutine for initializing general dense-block storage
    !-----------------------------------------------------------
    subroutine init_gen(self,idim,jdim,parent)
        class(denseblock_t), intent(inout)  :: self
        integer(ik),         intent(in)     :: idim, jdim, parent

        integer(ik) :: ierr

        ! Block parents
        self%parent_ = parent

        ! Allocate block storage
        ! Check if storage was already allocated and reallocate if necessary
        if (allocated(self%mat)) then
            deallocate(self%mat)
            allocate(self%mat(idim,jdim),stat=ierr)
!            if (ierr /= 0) call warn(3,__file__,__line__,'Allocation error')
!            print*, __line__
!            call warn(3,'Allocation error',__file__,__line__)
        else
            allocate(self%mat(idim,jdim),stat=ierr)
!            if (ierr /= 0) call warn(3,__file__,__line__,'Allocation error')
        end if

        ! Initialize to zero
        self%mat = 0._rk
    end subroutine



    !> Subroutine for initializing square dense-block storage
    !-----------------------------------------------------------
    subroutine init_square(self,bsize,parent)
        class(denseblock_t), intent(inout)  :: self
        integer(ik),         intent(in)     :: bsize, parent

        integer(ik) :: ierr

        ! Block parents
        self%parent_ = parent

        ! Allocate block storage
        ! Check if storage was already allocated and reallocate if necessary
        if (allocated(self%mat)) then
            deallocate(self%mat)
            allocate(self%mat(bsize,bsize),stat=ierr)
!            if (ierr /= 0) call warn(3,__file__,__line__,'Allocation error')
        else
            allocate(self%mat(bsize,bsize),stat=ierr)
!            if (ierr /= 0) call warn(3,__file__,__line__,'Allocation error')
        end if

        ! Initialize to zero
        self%mat = 0._rk
    end subroutine






    !> return i-dimension of block storage
    !------------------------------------------------------------
    function idim(self) result(i)
        class(denseblock_t), intent(in)   :: self
        integer(ik)                       :: i

        i = size(self%mat,1)
    end function

    !> return j-dimension of block storage
    !------------------------------------------------------------
    function jdim(self) result(j)
        class(denseblock_t), intent(in)   :: self
        integer(ik)                       :: j

        j = size(self%mat,2)
    end function

    !> return number of entries in block storage
    !------------------------------------------------------------
    function nentries(self) result(n)
        class(denseblock_t), intent(in)   :: self
        integer(ik)                  :: n

        n = size(self%mat,1) * size(self%mat,2)
    end function

    !> return index of block parent
    !------------------------------------------------------------
    function parent(self) result(par)
        class(denseblock_t), intent(in) :: self
        integer(ik)                     :: par

        par = self%parent_
    end function



    !> Resize dense-block storage
    !------------------------------------------------------------
    subroutine resize(self,idim,jdim)
        class(denseblock_t), intent(inout)  :: self
        integer(ik),         intent(in)     :: idim,jdim

        integer(ik) :: ierr

        ! Allocate block storage
        ! Check if storage was already allocated and reallocate if necessary
        if (allocated(self%mat)) then
            deallocate(self%mat)
            allocate(self%mat(idim,jdim),stat=ierr)
!            if (ierr /= 0) call warn(3,__file__,__line__,'Allocation error')
        else
            allocate(self%mat(idim,jdim),stat=ierr)
!            if (ierr /= 0) call warn(3,__file__,__line__,'Allocation error')
        end if

    end subroutine


    !> set index of parent
    !------------------------------------------------------------
    subroutine reparent(self,par)
        class(denseblock_t), intent(inout)  :: self
        integer(ik),         intent(in)     :: par

        self%parent_ = par
    end subroutine


    subroutine destructor(self)
        type(denseblock_t), intent(inout) :: self
        if (allocated(self%mat))    deallocate(self%mat)
    end subroutine

end module type_denseblock
