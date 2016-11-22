module precon_jacobi
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: DIAG
    use type_preconditioner,    only: preconditioner_t
    use type_chidgMatrix,       only: chidgMatrix_t
    use type_chidgVector,       only: chidgVector_t
    use type_chidg_data,        only: chidg_data_t
    use type_densematrix,       only: densematrix_t

    use mod_inv,                only: inv


    !> Block-Jacobi preconditioner
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!
    !-----------------------------------------------------------------------------
    type, extends(preconditioner_t) :: precon_jacobi_t

        !type(densematrix_t), allocatable    :: D(:,:,:)     !< inverse of block diagonal, (ndom,maxelems,ntime)
        type(chidgMatrix_t) :: D     !< inverse of block diagonal, (ndom,maxelems,ntime)

    contains
        procedure   :: init
        procedure   :: update
        procedure   :: apply

    end type precon_jacobi_t
    !******************************************************************************




contains


    !> Initialize preconditioner storage
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------
    subroutine init(self,data)
        class(precon_jacobi_t), intent(inout)   :: self
        type(chidg_data_t),     intent(in)      :: data

!        integer(ik) :: ndom, ntime
!        logical     :: increase_maxelems = .false.
!
!
!        ndom = data%ndomains()
!        ntime = data%ntime()
!
!
!        !
!        ! Get maximum number of elements
!        !
!        maxelems = 0
!        do idom = 1,ndom
!
!            increase_maxelems = ( data%mesh(idom)%nelem > maxelems )
!
!            if (increase_maxelems) then
!                maxelems = data%mesh(idom)%nelem
!            end if
!        
!        end do ! idom
!
!
!        !
!        ! Allocate storage
!        !
!        if (allocated(self%D)) deallocate(self%D)
!        allocate(self%D(ndom,maxelems,ntime), stat=ierr)
!        if (ierr /= 0) call AllocationError


        call self%D%init(data%mesh,mtype='Diagonal')


    end subroutine init
    !***************************************************************************************











    !> Compute the block diagonal inversion and store so it can be applied.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!
    !----------------------------------------------------------------------------------------
    subroutine update(self,A,b)
        class(precon_jacobi_t), intent(inout)   :: self
        type(chidgMatrix_t),    intent(in)      :: A
        type(chidgVector_t),    intent(in)      :: b


        integer(ik) :: idom, ielem, itime, diag


        !
        ! Copy the block diagonal from A to the preconditioner storage
        !
        do idom = 1,size(A%dom)
            do ielem = 1,size(A%dom(idom)%lblks,1)
                do itime = 1,size(A%dom(idom)%lblks,2)
                    diag = A%dom(idom)%lblks(ielem,itime)%get_diagonal()
                    self%D%dom(idom)%lblks(ielem,itime)%data_(1)%mat = A%dom(idom)%lblks(ielem,itime)%dmat(diag)
                    !self%D(idom,ielem,itime) = A%dom(idom)%lblks(ielem,DIAG)
                    !self%D(idom,ielem,itime) = A%dom(idom)%lblks(ielem,itime)%data_(index)
                end do
            end do
        end do



        !
        ! Replace the block diagonal D with Dinv
        !
        do idom = 1,size(A%dom)
            do ielem = 1,size(A%dom(idom)%lblks,1)
                do itime = 1,size(A%dom(idom)%lblks,2)
                    self%D%dom(idom)%lblks(ielem,itime)%data_(1)%mat = inv(self%D%dom(idom)%lblks(ielem,itime)%data_(1)%mat)
                end do
            end do
        end do


    end subroutine update
    !***************************************************************************************








    !> Apply the preconditioner to the krylov vector 'v' and return preconditioned vector 'z'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------
    function apply(self,A,v) result(z)
        class(precon_jacobi_t), intent(inout)   :: self
        type(chidgMatrix_t),    intent(in)      :: A
        type(chidgVector_t),    intent(in)      :: v

        type(chidgVector_t) :: z
        integer(ik)         :: idom, ielem, itime

        call self%timer%start()


        !
        ! Allocate 'z'
        !
        z = v
        call z%clear()


        !
        ! Apply Block-Diagonal preconditioner
        !
        do idom = 1,size(A%dom)
            do ielem = 1,size(A%dom(idom)%lblks,1)
                do itime = 1,size(A%dom(idom)%lblks,2)
                    z%dom(idom)%vecs(ielem)%vec = &
                        matmul(self%D%dom(idom)%lblks(ielem,itime)%data_(1)%mat,v%dom(idom)%vecs(ielem)%vec)
                end do
            end do
        end do


        call self%timer%stop()

    end function apply
    !****************************************************************************************













end module precon_jacobi
