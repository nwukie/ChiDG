module type_fgmres
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
    use operator_dot,       only: dot

    implicit none
        






    !> Generalized Minimum Residual linear system solver
    !!
    !!  @author Nathan A. Wukie
    !!
    !-------------------------------------------
    type, public, extends(matrixsolver_t) :: fgmres_t

        integer(ik) :: m = 100

    contains

        procedure   :: solve
    end type fgmres_t





contains


    !> Solution routine
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !--------------------------------------------------------------
    subroutine solve(self,A,x,b,M)
        class(fgmres_t),            intent(inout)               :: self
        type(blockmatrix_t),        intent(inout)               :: A
        type(blockvector_t),        intent(inout)               :: x
        type(blockvector_t),        intent(inout)               :: b
        class(preconditioner_t),    intent(inout), optional     :: M



        type(blockvector_t)                     :: r, r0, diff, xold, w, x0
        type(blockvector_t), allocatable        :: v(:), z(:)
        type(densematrix_t)                     :: D(size(x%lvecs))
        real(rk),            allocatable        :: h(:,:), h_square(:,:)
        real(rk),            allocatable        :: p(:), y(:), c(:), s(:), p_dim(:), y_dim(:)
        real(rk)                                :: pj, pjp, h_ij, h_ipj

        integer(ik) :: iparent, ierr, ivec, isol, nvecs, ielem
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
        ! Update preconditioner
        !
        call M%update(A,b)




        !
        ! Allocate and initialize Krylov vectors V
        !
        allocate(v(self%m+1),  &
                 z(self%m+1), stat=ierr)
        if (ierr /= 0) call AllocationError

        do ivec = 1,size(v)
            v(ivec) = b
            z(ivec) = b
            call v(ivec)%clear()
            call z(ivec)%clear()
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


        self%niter = 0


        res = 1._rk
        do while (res > self%tol)

            !
            ! Clear working variables
            !
            do ivec = 1,size(v)
                call v(ivec)%clear()
                call z(ivec)%clear()
            end do
            p = ZERO
            y = ZERO
            c = ZERO
            s = ZERO
            h = ZERO



            !
            ! Compute initial residual r0, residual norm, and normalized r0
            !
            r0      = self%residual(A,x0,b)
            v(1)    = r0/r0%norm()




            p(1) = r0%norm()
            !
            ! Outer GMRES Loop
            !
            nvecs = 0
            do j = 1,self%m
                nvecs = nvecs + 1
           


                !
                ! Apply preconditioner:  z(j) = Minv * v(j)
                !
                z(j) = M%apply(A,v(j))



                !
                ! Compute w = Av for the current iteration
                !
                w = A*z(j)



                !
                ! Orthogonalization loop
                !
                do i = 1,j

                    h(i,j) = dot(w,v(i))
                    
                    w  = w - h(i,j)*v(i)

                end do

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
                ! Update iteration counter
                !
                self%niter = self%niter + 1_ik



                !
                ! Test exit conditions
                !
                res = abs(p(j+1))
                print*, res
                converged = (res < self%tol)
                
                if ( converged ) then
                    exit
                end if




            end do  ! Outer GMRES Loop - m








            !
            ! Solve upper-triangular system y = hinv * p
            !
            if (allocated(h_square)) then
                deallocate(h_square,p_dim,y_dim)
            end if
            
            allocate(h_square(nvecs,nvecs), &
                    p_dim(nvecs),      &
                    y_dim(nvecs), stat=ierr)
            if (ierr /= 0) call AllocationError




            !
            ! Store h and p values to appropriately sized matrices
            !
            do l=1,nvecs
                do k=1,nvecs
                    h_square(k,l) = h(k,l)
                end do
                p_dim(l) = p(l)
            end do




            !
            ! Solve the system
            !
            h_square = inv(h_square)
            y_dim = matmul(h_square,p_dim)



            !
            ! Reconstruct solution
            !
            x = x0
            do isol = 1,nvecs
                x = x + y_dim(isol)*z(isol)
            end do



            !
            ! Test exit condition
            !
            if ( converged ) then
                exit
            else
                x0 = x
            end if



        end do   ! while









        !
        ! Report
        !
        err = self%error(A,x,b)
        print*, '   Matrix Solver Error: ', err

        call self%timer%stop()
        call self%timer%report('Matrix solver compute time: ')



    end subroutine solve



end module type_fgmres
