module mod_fgmres_standard
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: DIAG, ZERO, TWO
    use mod_inv,                only: inv
    implicit none
        









contains


    !> Solution routine
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !--------------------------------------------------------------
    subroutine fgmres_standard(A,x,b)
        real(rk), dimension(:,:),   intent(in)               :: A
        real(rk), dimension(:),     intent(inout)            :: x
        real(rk), dimension(:),     intent(in)               :: b



        real(rk), dimension(:), allocatable     :: r0, w, x0
        real(rk), allocatable                   :: v(:,:), z(:,:)
        real(rk),            allocatable        :: h(:,:), h_square(:,:)
        real(rk),            allocatable        :: p(:), y(:), c(:), s(:), p_dim(:), y_dim(:)
        real(rk)                                :: pj, pjp, h_ij, h_ipj

        integer(ik) :: ierr, isol, nvecs, mrestart
        integer(ik) :: i, j, k, l
        real(rk)    :: res, gam, tol

        logical     :: converged = .false.


        tol = 1.e-11_rk


        !
        ! Set GMRES parameters
        !
        mrestart = 200


        !
        ! Allocate and initialize Krylov vectors V
        !
        allocate(v(size(x),mrestart+1),  &
                 z(size(x),mrestart+1), stat=ierr)
        if (ierr /= 0) call AllocationError



        !
        ! Allocate hessenberg matrix to store orthogonalization
        !
        allocate(h(mrestart + 1, mrestart), stat=ierr)
        if (ierr /= 0) call AllocationError
        h = ZERO




        !
        ! Allocate vectors for solving hessenberg system
        !
        allocate(p(mrestart+1), &
                 y(mrestart+1), &
                 c(mrestart+1), &
                 s(mrestart+1), stat=ierr)
        if (ierr /= 0) call AllocationError
        p = ZERO
        y = ZERO
        c = ZERO
        s = ZERO



        !
        ! Set initial solution x. ZERO
        !
        x0 = x
        x0 = ZERO
        x  = ZERO




        res = 1._rk
        do while (res > tol)

            !
            ! Clear working variables
            !
            v = ZERO
            z = ZERO


            p = ZERO
            y = ZERO
            c = ZERO
            s = ZERO
            h = ZERO



            !
            ! Compute initial residual r0, residual norm, and normalized r0
            !
            r0      = b - matmul(A,x)

            
            v(:,1)    = r0/norm2(r0)




            p(1) = norm2(r0)
            !
            ! Outer GMRES Loop
            !
            nvecs = 0
            do j = 1, mrestart
                nvecs = nvecs + 1
           


                !
                ! Apply preconditioner:  z(j) = Minv * v(j)
                !
                !z(j) = M%apply(A,v(j))



                !
                ! Compute w = Av for the current iteration
                !
                w = matmul(A,v(:,j))



                !
                ! Orthogonalization loop
                !
                do i = 1,j

                    h(i,j) = dot_product(w,v(:,i))
                    
                    w  = w - h(i,j)*v(:,i)

                end do

                h(j+1,j) = norm2(w)




                !
                ! Compute next Krylov vector
                !
                v(:,j+1) = w/h(j+1,j)




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
                res = abs(p(j+1))
                !print*, '       ', res
                converged = (res < tol)
                
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
                !x = x + y_dim(isol)*z(isol)
                x = x + y_dim(isol)*v(:,isol)
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










    end subroutine fgmres_standard



end module mod_fgmres_standard
