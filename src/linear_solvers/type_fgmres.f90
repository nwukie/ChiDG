module type_fgmres
#include <messenger.h>
#include "petsc/finclude/petscvec.h"
    use petscvec,               only: VecAXPY, VecMDOT, VecMAXPY, VecDot, VecGetLocalVectorRead, VecRestoreLocalVectorRead, VecDuplicate, VecSet, VecDestroy

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
    use type_chidg_data,        only: chidg_data_t
    use type_chidg_vector
    use type_chidg_matrix

    use operator_chidg_dot,     only: dot
    use operator_chidg_mv,      only: chidg_mv, timer_comm, timer_blas
    implicit none


    !>  Flexible variant of the Generalized Minimum Residual algorithm for solving nonsymmetric
    !!  systems of linear equations.
    !!
    !!  Algorithms:
    !!      : Number of Krylov vectors       (&linear_solve  nkrylov=2000 /)
    !!      : Orthogonalization              (&linear_solve  orthogonalization='CGS' or 'MGS' /)
    !!      : Inner fgmres iteration   (&linear_solve  inner_fgmres=.true. /)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/11/2016
    !!
    !---------------------------------------------------------------------------------------------
    type, public, extends(linear_solver_t) :: fgmres_t

        ! Allocating storage on object so it can persist and the allocations can be reused
        type(chidg_vector_t)                                :: r, r0, w, x0, zr, deltaz
        type(chidg_vector_t),   allocatable, dimension(:)   :: v, z
        real(rk),               allocatable, dimension(:)   :: p, y, c, s, htmp
        real(rk),               allocatable, dimension(:,:) :: h

        type(timer_t)   :: timer_mv, timer_dot, timer_norm, timer_precon

        type(fgmres_t), pointer :: fgmres => null()

    contains
        procedure   :: solve
        procedure   :: allocate_krylov_storage
        procedure   :: print_report
        procedure   :: tear_down
        final       :: destroy_fgmres
    end type fgmres_t
    !*********************************************************************************************


contains





    !> Solution routine
    !!
    !!  NOTE: this subroutine is recursive so local variables should NOT be declared with the
    !!        SAVE attribute. The declaration and allocation of krylov vectors was relocated
    !!        to the object so that their allocation can be stored without using the SAVE 
    !!        attribute.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/11/2016
    !!
    !!  @author Nathan A. Wukie (AFRL) 
    !!  @date   6/23/2016
    !!  @note   parallelization
    !!
    !---------------------------------------------------------------------------------------------
    recursive subroutine solve(self,A,x,b,M,solver_controller,data)
        class(fgmres_t),            intent(inout)               :: self
        type(chidg_matrix_t),       intent(inout)               :: A
        type(chidg_vector_t),       intent(inout)               :: x
        type(chidg_vector_t),       intent(inout)               :: b
        class(preconditioner_t),    intent(inout), optional     :: M
        class(solver_controller_t), intent(inout), optional     :: solver_controller
        type(chidg_data_t),         intent(in),    optional     :: data

        integer(ik) :: ierr, ivec, nvecs, i, j
        real(rk)    :: res, r0norm, r0norm_initial, gam, pj, pjp, h_ij, h_ipj, norm_before, L_crit, crit
        logical     :: converged_absolute, converged_relative, converged_nmax, reorthogonalize
        real(rk),   allocatable :: y_dim(:)

        ! Inner fgmres iteration parameters
        if (self%inner_fgmres) then
            if (.not. associated(self%fgmres)) allocate(self%fgmres)
            self%fgmres%nkrylov           = self%inner_nkrylov
            self%fgmres%tol               = self%inner_tol
            self%fgmres%rtol              = self%inner_rtol
            self%fgmres%nmax              = self%inner_nmax
            self%fgmres%silence           = self%inner_silence
            self%fgmres%orthogonalization = self%inner_orthogonalization
            self%fgmres%inner_fgmres      = .false.
        end if


        ! Reset/Start timers
        call timer_comm%reset()
        call timer_blas%reset()
        call self%timer%reset()
        call self%timer%start()
        call write_line('           Linear Solver: ', io_proc=GLOBAL_MASTER, silence=(verbosity+self%silence<4))


        ! Update preconditioner
        if (present(solver_controller)) then
            if (solver_controller%update_preconditioner(A,M)) call M%update(A,b)
        else
            call M%update(A,b)
        end if


        ! Allocate storage
        call self%allocate_krylov_storage(b)
        associate(v=>self%v, z=>self%z, p=>self%p, y=>self%y, c=>self%c, s=>self%s, h=>self%h, &
                  htmp=>self%htmp, r=>self%r, r0=>self%r0, w=>self%w, zr=>self%zr, deltaz=>self%deltaz, x0=>self%x0)


        ! Set initial solution x. ZERO
        x0 = x
        deltaz = x
        call x0%clear()
        call x%clear()
        call deltaz%clear()


        self%niter = 0
        res = huge(1._rk)
        converged_absolute = .false.
        converged_relative = .false.
        converged_nmax     = .false.
        do while (res > self%tol)

            ! Clear working variables
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


            ! Compute initial residual r0, residual norm, and normalized r0
            r0     = self%residual(A,x0,b)
            r0norm = r0%norm(ChiDG_COMM)
            v(1)   = r0/r0norm
            p(1)   = r0norm


            ! Inner GMRES restart loop
            if (self%niter == 0) r0norm_initial = r0norm
            nvecs = 0
            do j = 1,self%nkrylov
                nvecs = nvecs + 1

                ! Apply preconditioner:  z(j) = Minv * v(j)
                call self%timer_precon%start()
                z(j) = M%apply(A,v(j))

                ! Inner fgmres correction
                if (self%inner_fgmres) then
                    zr = v(j) - chidg_mv(A,z(j))
                    call self%fgmres%solve(A,deltaz,zr,M,solver_controller)
                    z(j) = z(j) + deltaz
                end if
                call self%timer_precon%stop()

                ! Compute w = Av for the current iteration
                call self%timer_mv%start()
                w = chidg_mv(A,z(j))
                call self%timer_mv%stop()
                norm_before = w%norm(ChiDG_COMM)


                ! Orthogonalize once. Classical Gram-Schmidt
                call self%timer_dot%start()
                if (self%orthogonalization == 'CGS') then
                    call classical_gram_schmidt(w,v,j,h,htmp)
                else if (self%orthogonalization == 'MGS') then
                    call modified_gram_schmidt(w,v,j,h,htmp)
                end if
                call self%timer_dot%stop()
                ! End Orthogonalize once.

                ! Selective Reorthogonalization
                !
                ! Giraud and Langou
                ! "A robust criterion for the modified Gram-Schmidt algorithm with selective reorthogonalization."
                ! SIAM J. of Sci. Comp.     Vol. 25, No. 2, pp. 417-441.
                !
                ! They recommend L<1 for robustness, but it seems for these problems L can be increased.
                L_crit = 1.0_rk
                crit = sum(abs(h(1:j,j)))/norm_before
                reorthogonalize = (crit >= L_crit) 
                if (reorthogonalize) then
                    call self%timer_dot%start()
                    if (self%orthogonalization == 'CGS') then
                        call classical_gram_schmidt(w,v,j,h,htmp)
                    else if (self%orthogonalization == 'MGS') then
                        call modified_gram_schmidt(w,v,j,h,htmp)
                    end if
                    call self%timer_dot%stop()
                end if
                ! End Orthogonalize twice.


                ! Compute next Krylov vector
                call self%timer_norm%start()
                h(j+1,j) = w%norm(ChiDG_COMM)
                call self%timer_norm%stop()
                v(j+1) = w/h(j+1,j)


                ! Previous Givens rotations on h
                if (j /= 1) then
                    do i = 1,j-1
                        ! Need temp values here so we don't directly overwrite the h(i,j) and h(i+1,j) values 
                        h_ij     =  c(i)*h(i,j)  +  s(i)*h(i+1,j)
                        h_ipj    = -s(i)*h(i,j)  +  c(i)*h(i+1,j)

                        h(i,j)   = h_ij
                        h(i+1,j) = h_ipj
                    end do
                end if


                ! Compute next rotation
                gam  = hypot(h(j,j), h(j+1,j))
                c(j) = h(j,j)/gam
                s(j) = h(j+1,j)/gam


                ! Givens rotation on h
                h(j,j)   = gam
                h(j+1,j) = ZERO


                ! Givens rotation on p. Need temp values here so we aren't directly overwriting 
                ! the p(j) value until we want to
                pj  =  c(j)*p(j)
                pjp = -s(j)*p(j)

                p(j)   = pj
                p(j+1) = pjp


                ! Update iteration counter
                self%niter = self%niter + 1_ik


                ! Test exit conditions
                res = abs(p(j+1))
                call write_line(res, io_proc=GLOBAL_MASTER, silence=(verbosity+self%silence<4))

                converged_absolute = (res < self%tol)
                converged_relative = (res/r0norm_initial < self%rtol)
                converged_nmax     = ( self%niter >= self%nmax ) .and. (self%nmax > 0)
                
                ! Convergence check
                if ( converged_absolute .or. converged_relative .or. converged_nmax ) exit


            end do  ! Outer GMRES Loop - restarts after m iterations


            ! Solve upper-triangular system y = hinv * p
            h(1:nvecs,1:nvecs) = inv(h(1:nvecs,1:nvecs))
            y_dim = matmul(h(1:nvecs,1:nvecs),p(1:nvecs))


            ! Reconstruct solution
            x = x0
            do ivec = 1,nvecs
                x = x + y_dim(ivec)*z(ivec)
            end do

            ! Test exit condition
            if ( converged_absolute .or. converged_relative .or. converged_nmax ) exit

            ! Update initial solution for restart
            x0 = x

        end do ! while


        call self%print_report(A,x,b)

        end associate


    end subroutine solve
    !************************************************************************************************************




    !>  Classical Gram-Schmidt orthogonalization.
    !!
    !!  Given vector 'w' and set of vectors 'v(:)', procedure returns 'w' orthogonal
    !!  to vectors 'v(:)'. 
    !!
    !!  The Classical Gram-Schmidt algorithm achieves better parallel performance because
    !!  all dot products are syncronized using a single MPI_AllReduce, but the algorithm
    !!  is prone to numerical instability. Reorthogonalization is often used to prevent
    !!  breakdown due to inaccurate orthogonalization.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/11/2016
    !!  @date   12/27/2018 (AFRL) restructured to generalized orthogonalization routines.
    !!
    !----------------------------------------------------------------------------------------
    subroutine classical_gram_schmidt(w,v,j,h,htmp)
        type(chidg_vector_t),   intent(inout)   :: w
        type(chidg_vector_t),   intent(in)      :: v(:)
        integer(ik),            intent(in)      :: j
        real(rk),               intent(inout)   :: h(:,:)
        real(rk),               intent(inout)   :: htmp(:)

        integer(ik)    :: i, ierr
        real(rk)       :: dot_tmp(j)
        PetscErrorCode :: perr

        Vec :: local_vector_w, local_vector_v


        if (w%petsc_vector_created) then

            call chidg_signal(FATAL,'classical_gram_schmidt: not implemented for new petsc implementation.')

! This implementation does not seem to be working in parallel!
!            call VecDuplicate(w%petsc_vector, local_vector_w, perr)
!            if (perr /= 0) call chidg_signal(FATAL,'classical_gram_schmidt: error calling VecDuplicate.')
!            call VecDuplicate(v(1)%petsc_vector, local_vector_v, perr)
!            if (perr /= 0) call chidg_signal(FATAL,'classical_gram_schmidt: error calling VecDuplicate.')
!
!            call VecSet(local_vector_w, ZERO, perr)
!            call VecSet(local_vector_v, ZERO, perr)
!
!
!
!            call VecGetLocalVectorRead(w%petsc_vector, local_vector_w, perr)
!            if (perr /= 0) call chidg_signal(FATAL,'classical_gram_schmidt: error calling VecGetLocalVector.')
!
!            do i = 1,j
!                call VecGetLocalVectorRead(v(i)%petsc_vector, local_vector_v, perr)
!                if (perr /= 0) call chidg_signal(FATAL,'classical_gram_schmidt: error calling VecGetLocalVector.')
!                call VecDot(local_vector_w,local_vector_v,dot_tmp(i),perr)
!                if (perr /= 0) call chidg_signal(FATAL,'dot_comm: error calling petsc VecDot.')
!                call VecRestoreLocalVectorRead(v(i)%petsc_vector, local_vector_v, perr)
!                if (perr /= 0) call chidg_signal(FATAL,'classical_gram_schmidt: error calling VecRestoreLocalVector.')
!            end do
!
!
!            ! Reduce local dot-product values across processors, distribute result back to all
!            call MPI_AllReduce(dot_tmp,htmp,j,MPI_REAL8,MPI_SUM,ChiDG_COMM,ierr)
!
!            h(1:j,j) = h(1:j,j) + htmp(1:j)
!
!            ! Subtract off htmp(i)*v(i)
!            htmp = -htmp
!            do i = 1,j
!                !w = w - htmp(i)*v(i)
!                call VecAXPY(w%petsc_vector, htmp(i), v(i)%petsc_vector, perr)
!                if (perr /= 0) call chidg_signal(FATAL,'classical_gram_schmidt: error calling VecAXPY.')
!            end do
!
!            call VecRestoreLocalVectorRead(w%petsc_vector, local_vector_w, perr)
!            if (perr /= 0) call chidg_signal(FATAL,'classical_gram_schmidt: error calling VecRestoreLocalVector.')
!
!            call VecDestroy(local_vector_w, perr)
!            if (perr /= 0) call chidg_signal(FATAL,'classical_gram_schmidt: error calling VecDestroy.')
!            call VecDestroy(local_vector_v, perr)
!            if (perr /= 0) call chidg_signal(FATAL,'classical_gram_schmidt: error calling VecDestroy.')


! Previous implementation but doesn't work with pointers!
!            call VecMDot(w%petsc_vector, j, v(:)%petsc_vector, dot_tmp, perr)
!            if (perr /= 0) call chidg_signal(FATAL,'classical_gram_schmidt: error calling VecMDot.')
!
!            ! Add to h
!            h(1:j,j) = h(1:j,j) + dot_tmp(1:j)
!
!            ! Subtract from w
!            dot_tmp = -dot_tmp
!            call VecMAXPY(w%petsc_vector, j, dot_tmp, v(:)%petsc_vector, perr)
!            if (perr /= 0) call chidg_signal(FATAL,'classical_gram_schmidt: error calling VecMAXPY.')


        else

            do i = 1,j
                dot_tmp(i) = dot(w,v(i))
            end do

            ! Reduce local dot-product values across processors, distribute result back to all
            call MPI_AllReduce(dot_tmp,htmp,j,MPI_REAL8,MPI_SUM,ChiDG_COMM,ierr)

            h(1:j,j) = h(1:j,j) + htmp(1:j)

            do i = 1,j
                w = w - htmp(i)*v(i)
            end do

        end if

    end subroutine classical_gram_schmidt
    !****************************************************************************************



    !>  Modified Gram-Schmidt orthogonalization.
    !!
    !!  Given vector 'w' and set of vectors 'v(:)', procedure returns 'w' orthogonal
    !!  to vectors 'v(:)'. 
    !!
    !!  The Modified Gram-Schmidt algorithm performs poorly in parallel because
    !!  the algorithm syncronizes for each dot product. However, the Modified Gram-Schmidt
    !!  algorithm is more numerically stable than the Classical Gram-Schmidt algorithm.
    !!  Classical Gram-Schmidt with reorthogonalization is generally preferred for 
    !!  parallel efficiency reasons.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/11/2016
    !!  @date   12/27/2018 (AFRL) restructured to generalized orthogonalization routines.
    !!
    !----------------------------------------------------------------------------------------
    subroutine modified_gram_schmidt(w,v,j,h,htmp)
        type(chidg_vector_t),   intent(inout)   :: w
        type(chidg_vector_t),   intent(in)      :: v(:)
        integer(ik),            intent(in)      :: j
        real(rk),               intent(inout)   :: h(:,:)
        real(rk),               intent(inout)   :: htmp(:)

        integer(ik)     :: i
        PetscErrorCode  :: perr

        if (w%petsc_vector_created) then
            do i = 1,j
                htmp(i) = dot(w,v(i),ChiDG_COMM)
                call VecAXPY(w%petsc_vector,-htmp(i),v(i)%petsc_vector,perr)
            end do

        else 
            do i = 1,j
                htmp(i) = dot(w,v(i),ChiDG_COMM)
                w = w - htmp(i)*v(i)
            end do
        end if

        h(1:j,j) = h(1:j,j) + htmp(1:j)

    end subroutine modified_gram_schmidt
    !****************************************************************************************






    !>  Allocate krylov vector storage. Try to reuse allocation if possible to reduce
    !!  allocation overhead, since this might be used in inner-iteration.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   1/3/2019
    !!
    !----------------------------------------------------------------------------------------
    subroutine allocate_krylov_storage(self,b)
        class(fgmres_t),        intent(inout)   :: self
        type(chidg_vector_t),   intent(in)      :: b

        logical     :: allocate_vectors
        integer(ik) :: ivec, ierr

        ! Check if krylov storage needs allocated or reallocated 
        if (allocated(self%v)) then
            if ( size(self%v) /= (self%nkrylov+1) ) then
                call self%v%release()
                call self%z%release()
                deallocate(self%v,self%z,self%p,self%y,self%c,self%s,self%htmp,self%h)
                allocate_vectors = .true.
            end if
            allocate_vectors = .false.
        else
            allocate_vectors = .true.
        end if


        ! Allocate if triggered
        if (allocate_vectors) then
            allocate(self%v(self%nkrylov+1), &
                     self%z(self%nkrylov+1), &
                     self%p(self%nkrylov+1), &
                     self%y(self%nkrylov+1), &
                     self%c(self%nkrylov+1), &
                     self%s(self%nkrylov+1), &
                  self%htmp(self%nkrylov+1), &
                     self%h(self%nkrylov+1,self%nkrylov), stat=ierr)
            if (ierr /= 0) call AllocationError

            do ivec = 1,size(self%v)
                self%v(ivec) = b
                self%z(ivec) = b
            end do
        end if

    end subroutine allocate_krylov_storage
    !****************************************************************************************







    !>
    !!
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine print_report(self,A,x,b)
        class(fgmres_t),        intent(inout)   :: self
        type(chidg_matrix_t),   intent(inout)   :: A
        type(chidg_vector_t),   intent(inout)   :: x
        type(chidg_vector_t),   intent(in)      :: b

        real(rk) :: err

        err = self%error(A,x,b)
        call self%timer%stop()
        call write_line('   Linear Solver Error: ',         err,                  delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity+self%silence<4))
        call write_line('   Linear Solver compute time: ',  self%timer%elapsed(), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity+self%silence<4))
        call write_line('   Linear Solver Iterations: ',    self%niter,           delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity+self%silence<4))

        call write_line('   Preconditioner time: ',           self%timer_precon%elapsed(),                     delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity+self%silence<5))
        call write_line('       Precon time per iteration: ', self%timer_precon%elapsed()/real(self%niter,rk), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity+self%silence<5))

        call write_line('   MV time: ',                   self%timer_mv%elapsed(),                     delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity+self%silence<5))
        call write_line('       MV time per iteration: ', self%timer_mv%elapsed()/real(self%niter,rk), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity+self%silence<5))

        call write_line('   MV comm time: ',                   timer_comm%elapsed(),                     delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity+self%silence<5))
        call write_line('       MV comm time per iteration: ', timer_comm%elapsed()/real(self%niter,rk), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity+self%silence<5))
        call write_line('   MV blas time: ',                   timer_blas%elapsed(),                     delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity+self%silence<5))
        call write_line('       MV blas time per iteration: ', timer_blas%elapsed()/real(self%niter,rk), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity+self%silence<5))
        call timer_comm%reset()
        call timer_blas%reset()

        call write_line('   Dot time: ',  self%timer_dot%elapsed(),  delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity+self%silence<5))
        call write_line('   Norm time: ', self%timer_norm%elapsed(), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity+self%silence<5))

    end subroutine print_report
    !*******************************************************************************************




    !>
    !!
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------
    recursive subroutine tear_down(self)
        class(fgmres_t),    intent(inout)   :: self

        call self%r%release()
        call self%r0%release()
        call self%w%release()
        call self%x0%release()
        call self%zr%release()
        call self%deltaz%release()
        call self%v%release()
        call self%z%release()

        if (associated(self%fgmres)) then
            call self%fgmres%tear_down()
            deallocate(self%fgmres)
            self%fgmres => null()
        end if

    end subroutine tear_down
    !*******************************************************************************************










    !>
    !!
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------
    recursive subroutine destroy_fgmres(self)
        type(fgmres_t),    intent(inout)   :: self

        call self%tear_down()

    end subroutine destroy_fgmres
    !*******************************************************************************************


    

















end module type_fgmres
