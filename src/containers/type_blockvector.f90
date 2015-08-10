!> Data type for storing the matrix of dense blocks which hold the linearization for an algorithm
!!  @author Nathan A. Wukie
module type_blockvector
#include <messenger.h>
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: DIAG
    use type_mesh,          only: mesh_t
    use type_densevector,   only: densevector_t
    use DNAD_D
    implicit none


    type, public :: blockvector_t
        !> localblocks (nelem x 7)
        !!
        !!            densevector_t
        !!
        !!  elem #1:        vec #1
        !!  elem #2:        vec #2
        !!  elem #3:        vec #3
        !!    .
        !!    .
        type(densevector_t), allocatable :: lvecs(:)         !> Local element vectors
        integer(ik),         allocatable :: ldata(:,:)       !> Local block data     (nvars, nterms)


    contains
        !> Initializers
        generic,   public  :: init => init_vector   !> Initialize local vector
        procedure, private :: init_vector

        final :: destructor
    end type blockvector_t



    private
contains



    !> Subroutine for initializing local linearization matrix
    !-----------------------------------------------------------
    subroutine init_vector(self,mesh)
        class(blockvector_t), intent(inout)  :: self
        class(mesh_t),        intent(in)     :: mesh

        integer(ik) :: nelem, nblk, ierr, ielem, iblk, size1d, parent
        logical     :: new_elements

        nelem = mesh%nelem  !> Number of elements in the local block
        nblk  = 7           !> Number of blocks in the local linearization (1D => 3, 2D => 5, 3D => 7)

        ! ALLOCATE SIZE FOR 'lvecs'
        !----------------------------------------------------
        ! If vector was already allocated, deallocate and then reallocate vector size
        ! Reallocation would take place if the number of elements were changed
        if (allocated(self%lvecs)) then
            ! If the size is already allocated, check if the number of elements has changed.
            ! If so (new_elements), then reallocate matrix size.
            ! If not, do nothing
            new_elements = (mesh%nelem /= size(self%lvecs))
            if (new_elements) then
                deallocate(self%lvecs, self%ldata)
                allocate(self%lvecs(nelem), self%ldata(nelem,2), stat=ierr)
            end if
        else
            allocate(self%lvecs(nelem), self%ldata(nelem,2), stat=ierr)
        end if
        if (ierr /= 0) call AllocationError



        !> Loop through elements and call initialization for densevectors
        do ielem = 1,mesh%nelem
            size1d = mesh%elems(ielem)%neqns  *  mesh%elems(ielem)%nterms_s

            call self%lvecs(ielem)%init(size1d,parent)

            ! Store data about number of equations and number of terms in solution expansion
            self%ldata(ielem,1) = mesh%elems(ielem)%neqns
            self%ldata(ielem,2) = mesh%elems(ielem)%nterms_s
        end do



    end subroutine








    subroutine destructor(self)
        type(blockvector_t), intent(inout) :: self

    end subroutine

end module type_blockvector
