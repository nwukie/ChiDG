!> Data type for storing the matrix of dense blocks which hold the linearization for an algorithm
!!  @author Nathan A. Wukie
module type_blockmatrix
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: DIAG, ZERO
    use type_mesh,          only: mesh_t
    use type_densematrix,   only: densematrix_t
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
        type(densematrix_t), allocatable :: lblks(:,:)      !> Local domain blocks
        integer(ik),         allocatable :: ldata(:,:)      !> Local block data     (nvars, nterms)

!        type(denseblock_t), allocatable :: cblks(:,:)      !> Chimera inter-domain blocks


    contains
        !> Initializers
        generic,   public  :: init => initialize_linearization   !> Initialize local matrix
        procedure, private :: initialize_linearization


        !> Setters
        procedure :: store      !> Store linearization data
        procedure :: clear      !> Zero all data storage

        final :: destructor
    end type blockmatrix_t



    private
contains



    !> Subroutine for initializing local linearization matrix
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  mesh    mesh_t containing arrays of elements and faces
    !!
    !-----------------------------------------------------------
    subroutine initialize_linearization(self,mesh)
        class(blockmatrix_t), intent(inout)  :: self
        class(mesh_t),        intent(in)     :: mesh

        integer(ik) :: nelem, nblk, ierr, ielem, iblk, size1d, parent
        logical     :: new_elements

        nelem = mesh%nelem  !> Number of elements in the local block
        nblk  = 7           !> Number of blocks in the local linearization (1D => 3, 2D => 5, 3D => 7)

        !
        ! ALLOCATE SIZE FOR 'localblocks'
        ! If matrix was already allocated, deallocate and then reallocate matrix size
        ! Reallocation would take place if the number of elements were changed
        !
        if (allocated(self%lblks)) then
            ! If the size is already allocated, check if the number of elements has changed.
            ! If so (new_elements), then reallocate matrix size.
            ! If not, do nothing
            new_elements = (mesh%nelem /= size(self%lblks,1))
            if (new_elements) then
                deallocate(self%lblks, self%ldata)
                allocate(self%lblks(nelem,nblk), self%ldata(nelem,2), stat=ierr)
                if (ierr /= 0) call AllocationError
            end if
        else
            allocate(self%lblks(nelem,nblk), self%ldata(nelem,2), stat=ierr)
            if (ierr /= 0) call AllocationError
        end if



        !
        ! Loop through elements and call initialization for linearization denseblock matrices
        !
        do ielem = 1,mesh%nelem
            do iblk = 1,7
                size1d = mesh%elems(ielem)%neqns  *  mesh%elems(ielem)%nterms_s

                ! Parent is the element with respect to which the linearization is computed
                if (iblk == DIAG) then
                    parent = mesh%elems(ielem)%ielem
                else
                    parent = mesh%faces(ielem,iblk)%ineighbor
                end if

                ! Call initialization procedure if parent is not 0 (0 meaning there is no parent for that block, probably a boundary)
                if (parent /= 0) then
                    call self%lblks(ielem,iblk)%init(size1d,parent)

                    ! Store data about number of equations and number of terms in solution expansion
                    self%ldata(ielem,1) = mesh%elems(ielem)%neqns
                    self%ldata(ielem,2) = mesh%elems(ielem)%nterms_s
                end if
            end do
        end do

    end subroutine initialize_linearization








    !>  Stores derivative data to the linearization matrix
    !!
    !!      -- Given the integral data computed from the spatial discretization,
    !!         store the derivative values from the AD data types
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  integral    Array of modes from the spatial scheme, with embedded partial derivatives for the linearization matrix
    !!  @param[in]  ielem       Element for which the linearization was computed
    !!  @param[in]  iblk        Index of a block for the linearization of the given element
    !!  @param[in]  ivar        Index of the variable
    !!
    !------------------------------------------------------------------------------
    subroutine store(self,integral,ielem,iblk,ivar)
        class(blockmatrix_t),   intent(inout)   :: self
        type(AD_D),             intent(in)      :: integral(:)
        integer(ik),            intent(in)      :: ielem, iblk, ivar

        integer(ik) :: iarray, neqns, nterms, irow_start, irow


        ! Get stored information for the block
        neqns  = self%ldata(ielem,1)
        nterms = self%ldata(ielem,2)

        irow_start = ( (ivar - 1)  *  nterms)

        !> If sizes match, store derivative arrays to local block. Loop through integral values, for each value store its derivatives.
        !! The integral values here should be components of the RHS vector. An array of partial derivatives from an AD_D variable
        !! should be stored as a row in the block matrix.
        do iarray = 1,size(integral)
            !> Do a += operation to add derivatives to any that are currently stored
            irow = irow_start + iarray
            self%lblks(ielem,iblk)%mat(irow,:) = self%lblks(ielem,iblk)%mat(irow,:) + integral(iarray)%xp_ad_
        end do

    end subroutine store






    !>  Set all denseblock_t storage to zero
    !!
    !!  @author Nathan A. Wukie
    !--------------------------------------------------------------
    subroutine clear(self)
        class(blockmatrix_t),   intent(inout)   :: self

        integer(ik) :: ielem, iblk  ! do-loop counters

        ! For each element
        do ielem = 1,size(self%lblks,1)
            ! For each block linearization for the current element
            do iblk = 1,size(self%lblks,2)

                ! Check if the block storage is actually allocated
                if (allocated(self%lblks(ielem,iblk)%mat)) then
                    ! If so, set to ZERO
                    self%lblks(ielem,iblk)%mat = ZERO
                end if

            end do  ! iblk
        end do  ! ielem

    end subroutine clear






    !>  Builds full matrix representation of block matrix
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[inout] fullmat   Storage for full matrix representation
    !------------------------------------------------------------------
    subroutine build(self,fullmat)
        class(blockmatrix_t),       intent(inout)   :: self
        real(rk),   allocatable,    intent(inout)   :: fullmat(:,:)








    end subroutine build
























    subroutine destructor(self)
        type(blockmatrix_t), intent(inout) :: self

    end subroutine

end module type_blockmatrix
