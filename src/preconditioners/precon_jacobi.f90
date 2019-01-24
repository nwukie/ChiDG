module precon_jacobi
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: DIAG
    use type_preconditioner,    only: preconditioner_t
    use type_chidg_matrix,       only: chidg_matrix_t
    use type_chidg_vector,       only: chidg_vector_t
    use type_chidg_data,        only: chidg_data_t
    use type_densematrix,       only: densematrix_t

    use mod_inv,                only: inv

    use mod_fgmres_standard,    only: fgmres_standard


    !> Block-Jacobi preconditioner
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/24/2016
    !!
    !!
    !-----------------------------------------------------------------------------
    type, extends(preconditioner_t) :: precon_jacobi_t

        type(chidg_matrix_t) :: D     ! inverse of block diagonal, (ndom,maxelems,ntime)

    contains
        procedure   :: init
        procedure   :: update
        procedure   :: apply
        procedure   :: restrict

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
        type(chidg_matrix_t),    intent(in)      :: A
        type(chidg_vector_t),    intent(in)      :: b


        integer(ik) :: idom, ielem, itime, diag


        !
        ! Copy the block diagonal from A to the preconditioner storage
        !
        do idom = 1,size(A%dom)
            do ielem = 1,size(A%dom(idom)%lblks,1)
                do itime = 1,size(A%dom(idom)%lblks,2)
                    diag = A%dom(idom)%lblks(ielem,itime)%get_diagonal()
                    self%D%dom(idom)%lblks(ielem,itime)%data_(1)%mat = A%dom(idom)%lblks(ielem,itime)%dmat(diag)
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


        ! Update stamp
        call date_and_time(values=self%stamp)

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
    function apply(self,A,v,z_old) result(z)
        class(precon_jacobi_t), intent(inout)           :: self
        type(chidg_matrix_t),   intent(in)              :: A
        type(chidg_vector_t),   intent(in)              :: v
        type(chidg_vector_t),   intent(in), optional    :: z_old

        type(chidg_vector_t)    :: z
        integer(ik)             :: idom, ielem, itime, ndom
        real(rk),   allocatable :: mv(:)

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
                    !z%dom(idom)%vecs(ielem)%vec = &
                    !    matmul(self%D%dom(idom)%lblks(ielem,itime)%data_(1)%mat,v%dom(idom)%vecs(ielem)%vec)
                    mv = matmul(self%D%dom(idom)%lblks(ielem,itime)%data_(1)%mat,v%dom(idom)%vecs(ielem)%gettime(itime))
                    call z%dom(idom)%vecs(ielem)%settime(itime, mv)
                end do
            end do
        end do


        call self%timer%stop()

    end function apply
    !****************************************************************************************





    !>  Produce a restricted version of the current preconditioner.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   7/24/2017
    !!
    !!
    !-----------------------------------------------------------------------------------------
    function restrict(self,nterms_r) result(restricted)
        class(precon_jacobi_t), intent(in)  :: self
        integer(ik),            intent(in)  :: nterms_r

        type(precon_jacobi_t) :: restricted

        restricted%D = self%D%restrict(nterms_r)
        restricted%initialized = .true.

    end function restrict
    !****************************************************************************************











end module precon_jacobi
