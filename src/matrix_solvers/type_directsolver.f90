module type_directsolver
#include <messenger.h>
    use mod_kinds,          only: rk, ik
    use mod_inv,            only: inv
    use atype_matrixsolver, only: matrixsolver_t 
    use type_blockmatrix,   only: blockmatrix_t
    use type_blockvector,   only: blockvector_t        
        






    !> Direct solve of linear system of equations via 
    !! inversion of A in Ax=b
    !!
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------
    type, public, extends(matrixsolver_t) :: directsolver_t



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
    subroutine solve(self,A,x,b)
        class(directsolver_t),  intent(inout)   :: self
        type(blockmatrix_t),    intent(inout)   :: A
        type(blockvector_t),    intent(inout)   :: x
        type(blockvector_t),    intent(inout)   :: b



        real(rk),   allocatable                 :: Afull(:,:), Ainv(:,:)
        real(rk),   allocatable                 :: bfull(:), xfull(:)
        integer(ik) :: ierr


        !
        ! Build full-matrix/vector representation of A and b.
        !
        call A%build(Afull)

        allocate(Ainv(size(Afull,1),size(Afull,2)), &
                 bfull(size(Afull,1)),           & 
                 xfull(size(Afull,1)), stat=ierr)
        if (ierr /= 0) call AllocationError

        ! Build full-vector representation of b
        call b%build(bfull)


        !
        ! Invert A
        !
        Ainv = inv(Afull)



        !
        ! Compute x = Ainv * b
        !
        xfull = matmul(Ainv,bfull)


        !
        ! Distribute full-vector representation xfull to blockvector representation x
        !
        call x%distribute(xfull)


    end subroutine solve



end module
