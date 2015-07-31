!> Data type for storing the matrix of dense blocks which hold the linearization for an algorithm
!!  @author Nathan A. Wukie
module type_blockmatrix
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: DIAG
    use type_mesh,          only: mesh_t
    use type_denseblock,    only: denseblock_t
    use DNAD_D
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
        type(denseblock_t), allocatable :: lblks(:,:)       !> Local domain blocks
!        type(denseblock_t), allocatable :: cblks(:,:)      !> Chimera inter-domain blocks


    contains
        !> Initializers
        generic,   public  :: init => init_general, init_linearization   !> Initialize local matrix
        procedure, private :: init_general
        procedure, private :: init_linearization


        !> Setters
        procedure :: store      !> Store linearization data

        final :: destructor
    end type blockmatrix_t



    private
contains

    !> Subroutine for initializing local matrix size
    !-----------------------------------------------------------
    subroutine init_general(self,nrow,ncol)
        class(blockmatrix_t), intent(inout)  :: self
        integer(ik),          intent(in)     :: nrow,ncol

        integer(ik) :: ierr


        ! Allocate matrix size
        if (allocated(self%lblks)) then
            deallocate(self%lblks)
            allocate(self%lblks(nrow,ncol), stat=ierr)
            if (ierr /= 0) call AllocationError
        else
            allocate(self%lblks(nrow,ncol), stat=ierr)
            if (ierr /= 0) call AllocationError
        end if

    end subroutine


    !> Subroutine for initializing local linearization matrix
    !-----------------------------------------------------------
    subroutine init_linearization(self,mesh)
        class(blockmatrix_t), intent(inout)  :: self
        class(mesh_t),        intent(in)     :: mesh

        integer(ik) :: nelem, nblk, ierr, ielem, iblk, size1d, parent
        logical     :: new_elements

        nelem = mesh%nelem  !> Number of elements in the local block
        nblk  = 7           !> Number of blocks in the local linearization (1D => 3, 2D => 5, 3D => 7)

        ! ALLOCATE SIZE FOR 'localblocks'
        !----------------------------------------------------
        ! If matrix was already allocated, deallocate and then reallocate matrix size
        ! Reallocation would take place if the number of elements were changed
        if (allocated(self%lblks)) then
            ! If the size is already allocated, check if the number of elements has changed.
            ! If so (new_elements), then reallocate matrix size.
            ! If not, do nothing
            new_elements = (mesh%nelem /= size(self%lblks,1))
            if (new_elements) then
                deallocate(self%lblks)
                allocate(self%lblks(nelem,nblk), stat=ierr)
                if (ierr /= 0) call AllocationError
            end if
        else
            allocate(self%lblks(nelem,nblk), stat=ierr)
            if (ierr /= 0) call AllocationError
        end if




        !> Loop through elements and call initialization for linearization denseblock matrices
        do ielem = 1,mesh%nelem
            do iblk = 1,7
                size1d = mesh%elems(ielem)%neqns  *  mesh%elems(ielem)%nterms_s

                !> Parent is the element with which the linearization was computed
                if (iblk == DIAG) then
                    parent = mesh%elems(ielem)%ielem
                else
                    parent = mesh%faces(ielem,iblk)%ineighbor
                end if

                !> Call initialization procedure
                call self%lblks(ielem,iblk)%init(size1d,parent)
            end do
        end do



    end subroutine



    !>  Stores derivative data to the linearization matrix
    !!
    !!      -- Given the integral data computed from the spatial discretization,
    !!         store the derivative values from the AD data types
    !!
    !!  @author Nathan A. Wukie
    !!  @param[in]  integral    Array of modes from the spatial scheme, with embedded partial derivatives for the linearization matrix
    !!  @param[in]  ielem       Element for which the linearization was computed
    !!  @param[in]  iblk        Index of a block for the linearization of the given element
    !!  @param[in]  ivar        Index of the variable
    !!
    subroutine store(self,integral,ielem,iblk,ivar)
        class(blockmatrix_t),   intent(inout)   :: self
        type(AD_D),             intent(in)      :: integral(:)
        integer(ik),            intent(in)      :: ielem, iblk, ivar



    end subroutine





    subroutine destructor(self)
        type(blockmatrix_t), intent(inout) :: self

    end subroutine

end module type_blockmatrix
