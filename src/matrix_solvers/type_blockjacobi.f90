module type_blockjacobi
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_inv,                only: inv
    use atype_matrixsolver,     only: matrixsolver_t 
    use type_preconditioner,    only: preconditioner_t
    use type_blockmatrix,       only: blockmatrix_t
    use type_blockvector
        






    !> Direct solve of linear system of equations via 
    !! inversion of A in Ax=b
    !!
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------
    type, public, extends(matrixsolver_t) :: blockjacobi_t



    contains

        procedure   :: solve
    end type





contains


    !> Solution routine
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !--------------------------------------------------------------
    subroutine solve(self,A,x,b,M)
        class(blockjacobi_t),       intent(inout)           :: self
        type(blockmatrix_t),        intent(inout)           :: A
        type(blockvector_t),        intent(inout)           :: x
        type(blockvector_t),        intent(inout)           :: b
        class(preconditioner_t),    intent(inout), optional :: M



        type(blockvector_t)                     :: r, diff, xold
        integer(ik) :: iparent, ierr
        real(rk)    :: res

        integer(ik) :: i


        print*, '           Matrix Solver: '


        ! Replace the block diagonal D with inverse of D
        do ielem = 1,size(A%lblks,1)
            A%lblks(ielem,7)%mat = inv(A%lblks(ielem,7)%mat)
        end do


        xold = x
        call x%clear()
        call xold%clear()



        self%niter = 0
        res = 1._rk
        do while (res > self%tol)
            
            

            ! form b - Ax, except for diagonal
            r = b
            do ielem = 1,size(x%lvecs)
                do iblk = 1,6   ! don't multiply diagonal

                    if (allocated(A%lblks(ielem,iblk)%mat)) then
                        iparent = A%lblks(ielem,iblk)%parent()
                        r%lvecs(ielem)%vec = r%lvecs(ielem)%vec - matmul(A%lblks(ielem,iblk)%mat,xold%lvecs(iparent)%vec)
                    end if

                end do
            end do


            ! Multiply by inverse of block diagonal in A. This was already inverted 
            ! and restored to the block diagonal location in A.
            do ielem = 1,size(x%lvecs)
                x%lvecs(ielem)%vec = matmul(A%lblks(ielem,7)%mat,r%lvecs(ielem)%vec)
            end do



            diff = x - xold
            res = diff%norm()

            if (self%report) then 
                print*, res
            end if

            xold = x    ! update solution


            !
            ! Update iteration counter
            !
            self%niter = self%niter + 1

        end do ! res > self%tol












    end subroutine solve



end module type_blockjacobi
