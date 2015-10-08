module precon_jacobi
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: DIAG
    use type_preconditioner,    only: preconditioner_t
    use type_blockmatrix,       only: blockmatrix_t
    use type_blockvector,       only: blockvector_t
    use type_densematrix,       only: densematrix_t

    use mod_inv,    only: inv


    !> Block-Jacobi preconditioner
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------
    type, extends(preconditioner_t) :: precon_jacobi_t

        type(densematrix_t), allocatable    :: D(:)     !< inverse of block diagonal

    contains
        procedure   :: update
        procedure   :: apply

    end type




contains


    !> Compute the block diagonal inversion and store so it can be applied.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !---------------------------------------------------------------------
    subroutine update(self,A,b)
        class(precon_jacobi_t), intent(inout)   :: self
        type(blockmatrix_t),    intent(in)      :: A
        type(blockvector_t),    intent(in)      :: b


        integer(ik) :: ielem
        logical     :: reallocate   !< logical test for reallocating preconditioner storage

        
        !
        ! Allocate storage for the block-diagonal preconditioner
        !
        reallocate = ( size(A%lblks,1) /= size(self%D) ) 

        if (reallocate) then
            ! Deallocate if already allocated
            if (allocated(self%D)) deallocate(self%D)

            ! Allocate a denseblock for each element
            allocate(self%D(size(A%lblks,1)),stat=ierr)
            if (ierr /= 0) call AllocationError

        end if


        !
        ! Copy the block diagonal from A to the preconditioner storage
        !
        do ielem = 1,size(A%lblks,1)
            self%D(ielem) = A%lblks(ielem,DIAG)
        end do

        !
        ! Replace the block diagonal D with Dinv
        !
        do ielem = 1,size(A%lblks,1)
            self%D(ielem)%mat = inv(self%D(ielem)%mat)
        end do


    end subroutine update








    !> Apply the preconditioner to the krylov vector 'v' and return preconditioned vector 'z'
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !-------------------------------------------------------------------------
    function apply(self,A,v) result(z)
        class(precon_jacobi_t), intent(inout)   :: self
        type(blockmatrix_t),    intent(in)      :: A
        type(blockvector_t),    intent(in)      :: v

        type(blockvector_t) :: z
        integer(ik)         :: ielem


        !
        ! Allocate 'z'
        !
        z = v
        call z%clear()


        do ielem = 1,size(self%D)
            z%lvecs(ielem)%vec = matmul(self%D(ielem)%mat,v%lvecs(ielem)%vec)
        end do



    end function













end module
