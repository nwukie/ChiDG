module type_gmres
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: DIAG, ZERO, TWO
    use mod_inv,                only: inv
    use atype_matrixsolver,     only: matrixsolver_t 
    use type_preconditioner,    only: preconditioner_t
    use type_blockmatrix,       only: blockmatrix_t
    use type_densematrix,       only: densematrix_t
    use type_blockvector
    use operator_mv
    use operator_dot,           only: dot
    
    implicit none
        






    !> Generalized Minimum Residual linear system solver
    !!
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------
    type, public, extends(matrixsolver_t) :: gmres_t

        integer(ik) :: m = 1000


    contains

        procedure   :: solve
    end type gmres_t





contains


    !> Solution routine
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !--------------------------------------------------------------
    subroutine solve(self,A,x,b,M)
        class(gmres_t),             intent(inout)               :: self
        type(blockmatrix_t),        intent(inout)               :: A
        type(blockvector_t),        intent(inout)               :: x
        type(blockvector_t),        intent(inout)               :: b
        class(preconditioner_t),    intent(inout), optional     :: M



        type(blockvector_t)                     :: r, r0, diff, xold, w, x0
        type(blockvector_t), allocatable        :: v(:)
        real(rk),            allocatable        :: h(:,:), h_square(:,:)
        real(rk),            allocatable        :: p(:), y(:), c(:), s(:), p_dim(:), y_dim(:)
        real(rk)                                :: pj, pjp, h_ij, h_ipj

        integer(ik) :: iparent, ierr, ivec, isol, nvecs
        integer(ik) :: i, j, k, l                 ! Loop counters
        real(rk)    :: res, err, r0norm, gam

        logical     :: converged = .false.
        logical     :: max_iter  = .false.


        !
        ! Start timer
        !
        call self%timer%reset()
        call self%timer%start()
        print*, '           Matrix Solver: '




        !
        ! Allocate and initialize Krylov vectors V
        !
        allocate(v(self%m+1),stat=ierr)
        if (ierr /= 0) call AllocationError

        do ivec = 1,size(v)
            v(ivec) = b
            call v(ivec)%clear()
        end do




        !
        ! Allocate hessenberg matrix to store orthogonalization
        !
        allocate(h(self%m + 1, self%m), stat=ierr)
        if (ierr /= 0) call AllocationError
        h = ZERO




        !
        ! Allocate vectors for solving hessenberg system
        !
        allocate(p(self%m+1), &
                 y(self%m+1), &
                 c(self%m+1), &
                 s(self%m+1), stat=ierr)
        if (ierr /= 0) call AllocationError
        p = ZERO
        y = ZERO
        c = ZERO
        s = ZERO



        !
        ! Set initial solution x. ZERO
        !
        x0 = x
        call x0%clear()
        call x%clear()



        !
        ! Compute initial residual r0, residual norm, and normalized r0
        !
        r0      = self%residual(A,x,b)
        v(1)    = r0/r0%norm()




        p(1) = r0%norm()
        !
        ! Outer GMRES Loop
        !
        nvecs = 0
        do j = 1,self%m
            nvecs = nvecs + 1
       


            !
            ! Compute w = Av for the current iteration
            !
            w = A*v(j)


            !
            ! Orthogonalization loop
            !
            do i = 1,j

                h(i,j) = dot(w,v(i))
                
                w  = w - h(i,j)*v(i)


            end do  ! Inner GMRES Loop - i


            h(j+1,j) = w%norm()






            !
            ! Compute next Krylov vector
            !
            v(j+1) = w/h(j+1,j)










            !
            ! Previous Givens rotations on h
            !
            if (j /= 1) then
                do i = 1,j-1
                    ! Need temp values here so we don't directly overwrite the h(i,j) and h(i+1,j) values 
                    h_ij     = c(i)*h(i,j)  +  s(i)*h(i+1,j)
                    h_ipj    = -s(i)*h(i,j)  +  c(i)*h(i+1,j)


                    h(i,j)    = h_ij
                    h(i+1,j)  = h_ipj
                end do
            end if



            !
            ! Compute next rotation
            !
            gam  = sqrt( h(j,j)*h(j,j)  +  h(j+1,j)*h(j+1,j) )
            c(j) = h(j,j)/gam
            s(j) = h(j+1,j)/gam



            !
            ! Givens rotation on h
            !
            h(j,j)   = gam
            h(j+1,j) = ZERO


            !
            ! Givens rotation on p. Need temp values here so we aren't directly overwriting the p(j) value until we want to
            !
            pj = c(j)*p(j)
            pjp = -s(j)*p(j)


            p(j)     = pj
            p(j+1)   = pjp



            



            !
            ! Test exit conditions
            !
            print*, 'resid'
            print*, abs(p(j+1))
            converged = (abs(p(j+1)) < self%tol)
            
            if ( converged ) then
                exit
            end if




        end do  ! Outer GMRES Loop - j







        !
        ! Solve upper-triangular system y = hinv * p
        !
        allocate(h_square(nvecs,nvecs), &
                 p_dim(nvecs),      &
                 y_dim(nvecs), stat=ierr)
        if (ierr /= 0) call AllocationError




        ! Store h and p values to appropriately sized matrices
        print*, j
        do l=1,nvecs
            do k=1,nvecs
                h_square(k,l) = h(k,l)
            end do
            p_dim(l) = p(l)
        end do

       
        ! Solve the system
        h_square = inv(h_square)
        y_dim = matmul(h_square,p_dim)




        !
        ! Reconstruct solution
        !
        x = x0
        do isol = 1,nvecs
            x = x + y_dim(isol)*v(isol)
        end do





        err = self%error(A,x,b)
        print*, '   Matrix Solver Error: ', err




        !
        ! Report timings
        !
        call self%timer%stop()
        call self%timer%report('Matrix solver compute time: ')



    end subroutine solve



end module type_gmres
