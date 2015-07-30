!> Data type for storing the matrix of dense blocks which hold the linearization for an algorithm
!!  @author Nathan A. Wukie
module type_blockmatrix
    use messenger,          only: warn
    use mod_kinds,          only: rk,ik
    use type_denseblock,    only: denseblock_t
    implicit none


    type, public :: blockmatrix_t
        !> localblocks (nelem x 7)
        !!
        !!            xi_min   xi_max   eta_min   eta_max   zeta_min    zeta_max    diag
        !!
        !!  elem #1:
        !!  elem #2:
        !!  elem #3:
        !!    .
        !!    .
        type(denseblock_t), allocatable :: localblocks(:,:)
!        type(denseblock_t), allocatable :: chimerablocks(:,:)


    contains
        !> Initializers
        generic,   public  :: init_local => init_local_general, init_local_linearization   !> Initialize local matrix
        procedure, private :: init_local_general
        procedure, private :: init_local_linearization

        procedure :: nblocks   !> return number of matrix entries
        procedure :: idim       !> return i-dimension of matrix storage
        procedure :: jdim       !> return j-dimension of matrix storage

        !> Setters
        procedure :: resize     !> resize matrix storage

        final :: destructor
    end type blockmatrix_t



    private
contains

    !> Subroutine for initializing local matrix size
    !-----------------------------------------------------------
    subroutine init_local_general(self,nrow,ncol)
        class(blockmatrix_t), intent(inout)  :: self
        integer(ik),          intent(in)     :: nrow,ncol

        integer(ik) :: ierr


        ! Allocate matrix size
        if (allocated(self%localblocks) then
            deallocate(self%localblocks)
            allocate(self%localblocks(nrow,ncol,stat=ierr)
            if (ierr /= 0) call warn(3,__file__,__line__,'Allocation error')
        else
            allocate(self%localblocks(nrow,ncol,stat=ierr)
            if (ierr /= 0) call warn(3,__file__,__line__,'Allocation error')
        end if

    end subroutine


    !> Subroutine for initializing local linearization matrix
    !-----------------------------------------------------------
    subroutine init_local_linearization(self,mesh)
        class(blockmatrix_t), intent(inout)  :: self
        class(mesh_t),        intent(in)     :: mesh

        integer(ik) :: nelem,nblk,ierr
        logical     :: new_elements

        nelem = mesh%nelem  !> Number of elements in the local block
        nblk  = 7           !> Number of blocks in the local linearization (1D => 3, 2D => 5, 3D => 7)

        ! ALLOCATE SIZE FOR 'localblocks'
        !----------------------------------------------------
        ! If matrix was already allocated, deallocate and then reallocate matrix size
        ! Reallocation would take place if the number of elements were changed
        if (allocated(self%localblocks) then
            ! If the size is already allocated, check if the number of elements has changed.
            ! If so (new_elements), then reallocate matrix size.
            ! If not, do nothing
            new_elements = (mesh%nelem /= size(self%localblocks,1))
            if (new_elements) then
                deallocate(self%localblocks)
                allocate(self%localblocks(nelem,nblk,stat=ierr)
                if (ierr /= 0) call warn(3,__file__,__line__,'Allocation error')
            end if
        else
            allocate(self%localblocks(nelem,nblk,stat=ierr)
            if (ierr /= 0) call warn(3,__file__,__line__,'Allocation error')
        end if


        ! Loop through elements and initialize dense blocks
        do ielem_zeta
            do ielem_eta
                do ielem_xi




                end do  ! ielem_xi
            end do  ! ielem_eta
        end do  ! ielem_zeta




    end subroutine







    !> Subroutine for initializing dense-blocks
    !------------------------------------------------------------
    subroutine init_local_(self,mesh)
        class(blockmatrix_t),   intent(inout)   :: self)


    end subroutine



    !> return i-dimension of block storage
    !------------------------------------------------------------
    function idim(self) result(i)
        class(denseblock_t), intent(in)   :: self
        integer(ik)                       :: i

        i = size(self%vals,1)
    end function

    !> return j-dimension of block storage
    !------------------------------------------------------------
    function jdim(self) result(j)
        class(denseblock_t), intent(in)   :: self
        integer(ik)                       :: j

        j = size(self%vals,2)
    end function

    !> return number of entries in block storage
    !------------------------------------------------------------
    function nentries(self) result(n)
        class(denseblock_t), intent(in)   :: self
        integer(ik)                  :: n

        n = size(self%vals,1) * size(self%vals,2)
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
        if (allocated(self%vals)) then
            deallocate(self%vals)
            allocate(self%vals(idim,jdim),stat=ierr)
            if (ierr /= 0) call warn(3,__file__,__line__,'Allocation error')
        else
            allocate(self%vals(idim,jdim),stat=ierr)
            if (ierr /= 0) call warn(3,__file__,__line__,'Allocation error')
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
        deallocate(self%vals)
    end subroutine

end module type_blockmatrix
