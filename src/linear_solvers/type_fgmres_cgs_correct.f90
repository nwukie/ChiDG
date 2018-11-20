module type_fgmres_cgs_correct
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO
    use mod_inv,                only: inv
    use mod_chidg_mpi,          only: ChiDG_COMM, GLOBAL_MASTER
    use mod_io,                 only: verbosity
    use mpi_f08

    use type_timer,             only: timer_t
    use type_linear_solver,     only: linear_solver_t 
    use type_preconditioner,    only: preconditioner_t
    use type_solver_controller, only: solver_controller_t
    use type_fgmres_cgs,        only: fgmres_cgs_t
    use precon_jacobi,          only: precon_jacobi_t
    use type_chidg_data,        only: chidg_data_t
    use type_chidg_vector
    use type_chidg_matrix

    use operator_chidg_dot,     only: dot
    use operator_chidg_mv,      only: chidg_mv, timer_comm, timer_blas
    implicit none
        






    !> Generalized Minimum Residual linear system solver
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/11/2016
    !!
    !---------------------------------------------------------------------------------------------
    type, public, extends(linear_solver_t) :: fgmres_cgs_correct_t

        integer(ik) :: m = 2000
        !integer(ik) :: m = 10

    contains

        procedure   :: solve

    end type fgmres_cgs_correct_t
    !*********************************************************************************************





contains


    !> Solution routine
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/11/2016
    !!
    !!  @author Nathan A. Wukie (AFRL) 
    !!  @date   6/23/2016
    !!  @note   parallelization
    !!
    !---------------------------------------------------------------------------------------------
    subroutine solve(self,A,x,b,M,solver_controller,data)
        class(fgmres_cgs_correct_t),    intent(inout)               :: self
        type(chidg_matrix_t),           intent(inout)               :: A
        type(chidg_vector_t),           intent(inout)               :: x
        type(chidg_vector_t),           intent(inout)               :: b
        class(preconditioner_t),        intent(inout), optional     :: M
        class(solver_controller_t),     intent(inout), optional     :: solver_controller
        type(chidg_data_t),             intent(in),    optional     :: data

        type(timer_t)   :: timer_mv, timer_dot, timer_norm, timer_precon


        type(chidg_vector_t)                :: r, r0, diff, xold, w, x0, zr, deltaz
        type(chidg_vector_t),   allocatable :: v(:), z(:)
        real(rk),               allocatable :: h(:,:), h_square(:,:), dot_tmp(:), htmp(:,:)
        real(rk),               allocatable :: p(:), y(:), c(:), s(:), p_dim(:), y_dim(:)
        real(rk)                            :: pj, pjp, h_ij, h_ipj, norm_before, norm_after, L_crit, crit

        integer(ik) :: iparent, ierr, ivec, isol, nvecs, ielem, xstart, xend, icorrect, idom, itime, diag, ismooth
        integer(ik) :: i, j, k, l, ii, ih                 ! Loop counters
        real(rk)    :: res, err, r0norm, gam, delta 

        logical :: converged = .false.
        logical :: max_iter  = .false.
        logical :: reorthogonalize = .false.



        type(fgmres_cgs_t) :: linear_solver
!        type(precon_jacobi_t)   :: jacobi

        linear_solver%m = 100
        !linear_solver%tol = 1.e-1_rk
        !linear_solver%rtol = 6.e-1_rk
        linear_solver%tol = 1.e-1_rk
        linear_solver%rtol = 1.e-1_rk
        linear_solver%maxiter = 100

        
        !
        !   BEGIN JACOBI Smoother
        !
!        call jacobi%init(data)
!        call jacobi%update(A,b)





        !
        ! Reset/Start timers
        !
        call timer_comm%reset()
        call timer_blas%reset()

        call self%timer%reset()
        call self%timer%start()
        call write_line('           Linear Solver: ', io_proc=GLOBAL_MASTER, silence=(verbosity<4))



        !
        ! Update preconditioner
        !
        if (present(solver_controller)) then
            if (solver_controller%update_preconditioner(A)) call M%update(A,b)
        else
            call M%update(A,b)
        end if


        

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
        allocate(h(self%m + 1, self%m), dot_tmp(self%m+1), htmp(self%m + 1, self%m), stat=ierr)
        if (ierr /= 0) call AllocationError
        h       = ZERO
        htmp    = ZERO
        dot_tmp = ZERO



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
        zr = x
        deltaz = x
        call x0%clear()
        call x%clear()
        call zr%clear()
        call deltaz%clear()


        self%niter = 0



        !res = 1000000000000._rk
        res = huge(1._rk)
        do while (res > self%tol)

            !
            ! Clear working variables
            !
            do ivec = 1,size(v)
                call v(ivec)%clear()
                call z(ivec)%clear()
            end do


            p    = ZERO
            y    = ZERO
            c    = ZERO
            s    = ZERO
            h    = ZERO
            htmp = ZERO



            !
            ! Compute initial residual r0, residual norm, and normalized r0
            !
            r0     = self%residual(A,x0,b)
            r0norm = r0%norm(ChiDG_COMM)
            v(1)   = r0/r0norm
            p(1)   = r0norm


            !
            ! Inner GMRES restart loop
            !
            nvecs = 0
            do j = 1,self%m
                nvecs = nvecs + 1


                !
                ! Apply preconditioner:  z(j) = Minv * v(j)
                !
                call timer_precon%start()
                z(j) = M%apply(A,v(j))
                call timer_precon%stop()


                !! Compute residual and use GMRES inner iterations to compute 
                !! approximate correction
                zr = v(j) - chidg_mv(A,z(j))
                call linear_solver%solve(A,deltaz,zr,M,solver_controller)
                z(j) = z(j) + deltaz

!                do ismooth = 1,3
!                    zr = v(j) - chidg_mv(A,z(j))
!                    deltaz = jacobi%apply(A,zr)
!                    z(j) = z(j) + 0.66_rk*deltaz
!                end do







                !
                ! Compute w = Av for the current iteration
                !
                call timer_mv%start()
                w = chidg_mv(A,z(j))
                call timer_mv%stop()



                norm_before = w%norm(ChiDG_COMM)



                call timer_dot%start()
                !
                ! Orthogonalize once. Classical Gram-Schmidt
                !
                do i = 1,j
                    ! Compute the local dot product
                    dot_tmp(i) = dot(w,v(i))
                end do

                ! Reduce local dot-product values across processors, distribute result back to all
                !call MPI_AllReduce(dot_tmp,h(:,j),j,MPI_REAL8,MPI_SUM,ChiDG_COMM,ierr)
                call MPI_AllReduce(dot_tmp,h(1:j,j),j,MPI_REAL8,MPI_SUM,ChiDG_COMM,ierr)
                

                do i = 1,j
                    w = w - h(i,j)*v(i)
                end do
                call timer_dot%stop()




                call timer_norm%start()
                h(j+1,j) = w%norm(ChiDG_COMM)
                norm_after = h(j+1,j)
                call timer_norm%stop()
                !
                ! End Orthogonalize once.
                !



                !
                ! Selective reorthogonalization
                !
                ! Giraud and Langou
                ! "A robust criterion for the modified Gram-Schmidt algorithm with selective reorthogonalization."
                ! SIAM J. of Sci. Comp.     Vol. 25, No. 2, pp. 417-441.
                !
                ! They recommend L<1 for robustness, but it seems for these problems L can be increased.
                !
                L_crit = 1.0_rk
                crit = sum(abs(h(1:j,j)))/norm_before

                !
                ! Orthogonalize twice. Classical Gram-Schmidt
                !
                reorthogonalize = (crit >= L_crit) 
                if (reorthogonalize) then
                    do i = 1,j
                        ! Compute the local dot product
                        dot_tmp(i) = dot(w,v(i))
                    end do


                    ! Reduce local dot-product values across processors, distribute result back to all
                    !call MPI_AllReduce(dot_tmp,htmp(:,j),j,MPI_REAL8,MPI_SUM,ChiDG_COMM,ierr)
                    call MPI_AllReduce(dot_tmp,htmp(1:j,j),j,MPI_REAL8,MPI_SUM,ChiDG_COMM,ierr)
                    !h(:,j) = h(:,j) + htmp(:,j)


                    !htmp(1,1) = 5.0

                    h = htmp + h
                    !do ih = 1,size(h,2)
                    !    h(:,ih) = h(:,ih) + htmp(:,ih)
                    !end do


                    do i = 1,j
                        w = w - htmp(i,j)*v(i)
                    end do



                    call timer_norm%start()
                    h(j+1,j) = w%norm(ChiDG_COMM)
                    call timer_norm%stop()
                end if
                !
                ! End Orthogonalize twice.
                !
                



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
                        h_ij     =  c(i)*h(i,j)  +  s(i)*h(i+1,j)
                        h_ipj    = -s(i)*h(i,j)  +  c(i)*h(i+1,j)


                        h(i,j)   = h_ij
                        h(i+1,j) = h_ipj
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
                pj  =  c(j)*p(j)
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
                call write_line(res, io_proc=GLOBAL_MASTER, silence=(verbosity<4))
                converged = (res < self%tol)
                
                if ( converged ) then
                    exit
                end if


            end do  ! Outer GMRES Loop - restarts after m iterations





            !
            ! Solve upper-triangular system y = hinv * p
            !
            if (allocated(h_square)) then
                deallocate(h_square,p_dim,y_dim)
            end if
            
            allocate(h_square(nvecs,nvecs), &
                     p_dim(nvecs),          &
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
        call self%timer%stop()
        call write_line('   Linear Solver Error: ',         err,                  delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<4))
        call write_line('   Linear Solver compute time: ',  self%timer%elapsed(), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<4))
        call write_line('   Linear Solver Iterations: ',    self%niter,           delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<4))

        !call self%timer%report('Linear solver compute time: ')
        !call timer_precon%report('Preconditioner time: ')
        call write_line('   Preconditioner time: ',           timer_precon%elapsed(),                     delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))
        call write_line('       Precon time per iteration: ', timer_precon%elapsed()/real(self%niter,rk), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))


        !call timer_mv%report('MV time: ')
        call write_line('   MV time: ',                   timer_mv%elapsed(),                     delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))
        call write_line('       MV time per iteration: ', timer_mv%elapsed()/real(self%niter,rk), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))

        !call timer_comm%report('MV comm time: ')
        call write_line('   MV comm time: ',                   timer_comm%elapsed(),                     delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))
        call write_line('       MV comm time per iteration: ', timer_comm%elapsed()/real(self%niter,rk), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))
        !call timer_blas%report('MV blas time: ')
        call write_line('   MV blas time: ',                   timer_blas%elapsed(),                     delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))
        call write_line('       MV blas time per iteration: ', timer_blas%elapsed()/real(self%niter,rk), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))
        call timer_comm%reset()
        call timer_blas%reset()

        !call timer_dot%report('Dot time: ')
        !call timer_norm%report('Norm time: ')
        call write_line('   Dot time: ',  timer_dot%elapsed(),  delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))
        call write_line('   Norm time: ', timer_norm%elapsed(), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))



    end subroutine solve
    !************************************************************************************************************



end module type_fgmres_cgs_correct
