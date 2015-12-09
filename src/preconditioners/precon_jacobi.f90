module precon_jacobi
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: DIAG
    use type_preconditioner,    only: preconditioner_t
    use type_chidgMatrix,       only: chidgMatrix_t
    use type_chidgVector,       only: chidgVector_t
    use type_chidg_data,        only: chidg_data_t
    use type_densematrix,       only: densematrix_t

    use mod_inv,    only: inv
    !use mod_gaussseidel_standard,   only: gaussseidel_standard
    use mod_fgmres_standard,   only: fgmres_standard


    !> Block-Jacobi preconditioner
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------
    type, extends(preconditioner_t) :: precon_jacobi_t

        type(densematrix_t), allocatable    :: D(:,:)     !< inverse of block diagonal, (ndom,maxelems)

    contains
        procedure   :: init
        procedure   :: update
        procedure   :: apply

    end type




contains


    !> Initialize preconditioner storage
    !!
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------
    subroutine init(self,data)
        class(precon_jacobi_t), intent(inout)   :: self
        type(chidg_data_t),     intent(in)      :: data

        integer(ik) :: ndom
        logical     :: increase_maxelems = .false.


        ndom = data%ndomains()


        !
        ! Get maximum number of elements
        !
        maxelems = 0
        do idom = 1,ndom

            increase_maxelems = ( data%mesh(idom)%nelem > maxelems )

            if (increase_maxelems) then
                maxelems = data%mesh(idom)%nelem
            end if
        
        end do ! idom


        !
        ! Allocate storage
        !
        if (allocated(self%D)) deallocate(self%D)
        allocate(self%D(ndom,maxelems), stat=ierr)
        if (ierr /= 0) call AllocationError






    end subroutine init











    !> Compute the block diagonal inversion and store so it can be applied.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !---------------------------------------------------------------------
    subroutine update(self,A,b)
        class(precon_jacobi_t), intent(inout)   :: self
        type(chidgMatrix_t),    intent(in)      :: A
        type(chidgVector_t),    intent(in)      :: b


        integer(ik) :: ielem, ndom
        logical     :: reallocate   !< logical test for reallocating preconditioner storage
        type(densematrix_t) :: junk

        ndom = size(A%dom)

        !
        ! Copy the block diagonal from A to the preconditioner storage
        !
        do idom = 1,ndom
            do ielem = 1,size(A%dom(idom)%lblks,1)
                self%D(idom,ielem) = A%dom(idom)%lblks(ielem,DIAG)
            end do
        end do

        !
        ! Junk call to compute condition numbers
        !
        !do idom = 1,ndom
        !    do ielem = 1,size(A%dom(idom)%lblks,1)
        !        print*, 'Invert, display condition number for element: ', ielem
        !        junk%mat = inv(A%dom(idom)%lblks(ielem,DIAG)%mat)
        !    end do
        !end do




        !
        ! Replace the block diagonal D with Dinv
        !
        !do idom = 1,ndom
        !    do ielem = 1,size(A%dom(idom)%lblks,1)
        !        self%D(idom,ielem)%mat = inv(self%D(idom,ielem)%mat)
        !    end do
        !end do


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
        type(chidgMatrix_t),    intent(in)      :: A
        type(chidgVector_t),    intent(in)      :: v

        type(chidgVector_t) :: z
        integer(ik)         :: ielem, idom, ndom


        ndom = size(A%dom)

        !
        ! Allocate 'z'
        !
        z = v
        call z%clear()


        !
        ! Apply Block-Diagonal preconditioner
        !
        do idom = 1,ndom
            do ielem = 1,size(A%dom(idom)%lblks,1)
                call fgmres_standard(self%D(idom,ielem)%mat, z%dom(idom)%lvecs(ielem)%vec, v%dom(idom)%lvecs(ielem)%vec)
                !z%dom(idom)%lvecs(ielem)%vec = matmul(self%D(idom,ielem)%mat,v%dom(idom)%lvecs(ielem)%vec)
            end do
        end do



    end function













end module
