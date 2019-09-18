module type_jfnk
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, DIAG
    use mod_hdfio,              only: write_fields_hdf
    use mod_chidg_mpi,          only: ChiDG_COMM, GLOBAL_MASTER, IRANK, NRANK
    use mod_io,                 only: verbosity
    use mod_inv,                only: inv
    use mpi_f08,                only: MPI_Barrier

    use type_chidg_data,        only: chidg_data_t
    use type_nonlinear_solver,  only: nonlinear_solver_t
    use type_system_assembler,  only: system_assembler_t
    use type_linear_solver,     only: linear_solver_t
    use type_preconditioner,    only: preconditioner_t
    use type_solver_controller, only: solver_controller_t
    use type_chidg_vector

    use mod_modify_jacobian_sst, only: modify_jacobian_sst
    use ieee_arithmetic,        only: ieee_is_nan
    use operator_chidg_dot,     only: dot
    use operator_chidg_mv,      only: chidg_mv, timer_comm, timer_blas
    implicit none
    private



    !>  Solution advancement via the backward-euler method
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !------------------------------------------------------------------------------------------
    type, extends(nonlinear_solver_t), public :: jfnk_t


    contains
    
        procedure   :: solve
        !procedure   :: linear_residual 
        final       :: destructor

    end type jfnk_t
    !******************************************************************************************










contains





    !>  Solve for update 'dq'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine solve(self,data,system,linear_solver,preconditioner,solver_controller)
        class(jfnk_t),                          intent(inout)           :: self
        type(chidg_data_t),                     intent(inout)           :: data
        class(system_assembler_t),  optional,   intent(inout)           :: system
        class(linear_solver_t),     optional,   intent(inout)           :: linear_solver
        class(preconditioner_t),    optional,   intent(inout),  target  :: preconditioner
        class(solver_controller_t), optional,   intent(inout),  target  :: solver_controller

        character(100)          :: filename
        integer(ik)             :: itime, nsteps, ielem, iblk, iindex,  &
                                   niter, ieqn, idom, ierr,                     &
                                   rstart, rend, cstart, cend, nterms, imat, iwrite, step, eqn_ID, icfl

        real(rk)                :: dtau, amp, cfl, timing, resid, resid_prev, resid0, resid_new,    &
                                   alpha, f0, fn, forcing_term, residual_ratio
        real(rk), allocatable   :: vals(:), cfln(:), rnorm0(:), rnorm(:)
        type(chidg_vector_t)    :: b, qn, qold, qnew, dqdtau, q0, x0, w, r0
        logical                 :: searching, absolute_convergence, relative_convergence, update_preconditioner

        type(solver_controller_t),  target  :: default_controller
        class(solver_controller_t), pointer :: controller

        type(chidg_vector_t),   allocatable :: v(:), z(:)
        real(rk),               allocatable :: h(:,:), h_square(:,:), dot_tmp(:), htmp(:,:)
        real(rk),               allocatable :: p(:), y(:), c(:), s(:), p_dim(:), y_dim(:)
        real(rk)                            :: pj, pjp, h_ij, h_ipj, norm_before, norm_after, L_crit, crit

        integer(ik) :: iparent, ivec, isol, nvecs, xstart, xend
        integer(ik) :: gmres_niter, gmres_m, i, j, k, l, ii, ih                 ! Loop counters
        real(rk)    :: res, err, r0norm, gam, delta, eps

        logical :: converged = .false.
        logical :: max_iter  = .false.
        logical :: reorthogonalize = .false.



        if (present(solver_controller)) then
            controller => solver_controller
        else
            controller => default_controller
        end if


        gmres_m = 2000
      
        call write_line('NONLINEAR SOLVER', io_proc=GLOBAL_MASTER, silence=(verbosity<2))
        call write_line("iter","|R(Q)|","CFL", "Linear Solver(niter)", "LHS Updated", "Preconditioner Updated", delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))


        associate ( q   => data%sdata%q,    &
                    dq  => data%sdata%dq,   &
                    rhs => data%sdata%rhs,  &
                    lhs => data%sdata%lhs)

        allocate(v(gmres_m+1),  &
                 z(gmres_m+1), stat=ierr)
        if (ierr /= 0) call AllocationError

        do ivec = 1,size(v)
            v(ivec) = q 
            z(ivec) = q 
            call v(ivec)%clear()
            call z(ivec)%clear()
        end do

        allocate(h(gmres_m + 1, gmres_m), dot_tmp(gmres_m+1), htmp(gmres_m + 1, gmres_m), stat=ierr)
        if (ierr /= 0) call AllocationError

        allocate(p(gmres_m+1), &
                 y(gmres_m+1), &
                 c(gmres_m+1), &
                 s(gmres_m+1), stat=ierr)
        if (ierr /= 0) call AllocationError

            !
            ! start timer
            !
            !call self%timer%reset()
            !call self%timer%start()


            !
            ! Startup values
            !
            absolute_convergence = .true.
            relative_convergence = .true.
            update_preconditioner = .true.
            qn     = q      ! Store qn, since q will be operated on in the inner loop
            resid  = ONE    ! Force inner loop entry
            niter  = 0      ! Initialize inner loop counter


            !
            ! NONLINEAR CONVERGENCE INNER LOOP
            !
            do while ( absolute_convergence .and. relative_convergence )
                niter = niter + 1

                if (niter==1) then
                    update_preconditioner=.true.
                elseif (modulo(niter,10)==0) then
                    update_preconditioner=.true.
                else
                    update_preconditioner=.false.
                end if


                !
                ! Store the value of the current inner iteration solution (k) 
                ! for the solution update (n+1), q_(n+1)_k
                !
                qold = q


                !
                ! Update Spatial Residual and Linearization (rhs, lin)
                !
                !if (present(solver_controller)) then
                !    if ( niter <= 2)  then
                !        residual_ratio = ONE
                !    else
                !        residual_ratio = resid/resid_prev
                !    end if

                !    call system%assemble( data,             &
                !                          timing=timing,    &
                !                          differentiate=solver_controller%update_lhs(niter,residual_ratio) )
                !else
                !    call system%assemble(data,timing=timing,differentiate=.true.)
                !end if
                if ( niter <= 2)  then
                    residual_ratio = ONE
                else
                    residual_ratio = resid/resid_prev
                end if

                call system%assemble( data,             &
                                      timing=timing,    &
                                      differentiate=update_preconditioner)!controller%update_lhs(lhs,niter,residual_ratio) )
                resid_prev = resid
                resid      = rhs%norm(ChiDG_COMM)



                !
                ! Print diagnostics, check tolerance.
                !
                call self%residual_norm%push_back(resid)
                call write_line('|R| = ', resid, io_proc=GLOBAL_MASTER, silence=(verbosity<4))
                if ( resid < self%tol ) then
                    call write_line(niter, resid, cfln(1), 0, delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))
                    exit
                end if


                if ( ieee_is_nan(resid) ) then
                    call chidg_signal(FATAL,"jfnk%solve: NaN residual calculation. Check initial solution and operator objects.")
                end if

                call self%residual_time%push_back(timing)


                !
                ! Compute and store residual norm for each field
                !
                if (niter == 1) then
                    resid0 = rhs%norm(ChiDG_COMM)
                    rnorm0 = rhs%norm_fields(ChiDG_COMM)
                    rnorm0 = resid0 !override
                end if

                rnorm = rhs%norm_fields(ChiDG_COMM)
                rnorm = resid !override




                !
                ! During the first few iterations, allow the initial residual norm to update
                ! if it has increased. Otherwise, if a solution converged to an error floor
                ! and the residual raised a little bit, then the CFL would essentially reset
                ! from infinity to CFL0, which we do not want.
                !
                if (niter < 5) then
                    where (rnorm > rnorm0)
                        rnorm0 = rnorm0 + 0.3_rk*(rnorm - rnorm0)
                    end where
                end if



                !
                ! Compute new cfl for each field
                !
                if (allocated(cfln)) deallocate(cfln)
                allocate(cfln(size(rnorm)), stat=ierr)
                if (ierr /= 0) call AllocationError

                where (rnorm /= 0.)
                    cfln = self%cfl0*(rnorm0/rnorm)
                else where
                    cfln = 0.1
                end where

                !if (IRANK == GLOBAL_MASTER) then
                !    call add_to_line("  CFL: ")
                !    do icfl = 1,size(cfln)
                !        call add_to_line(cfln(icfl))
                !    end do
                !    call send_line()
                !end if



                !
                ! Compute element-local pseudo-timestep
                !
                call compute_pseudo_timestep(data,cfln)


                !
                ! Add mass/dt to sub-block diagonal in dR/dQ
                !
                if (update_preconditioner) then
                do idom = 1,data%mesh%ndomains()
                    do ielem = 1,data%mesh%domain(idom)%nelem
                        eqn_ID = data%mesh%domain(idom)%elems(ielem)%eqn_ID
                        do itime = 1,data%mesh%domain(idom)%ntime
                            do ieqn = 1,data%eqnset(eqn_ID)%prop%nprimary_fields()

                                ! get element-local timestep
                                dtau = data%mesh%domain(idom)%elems(ielem)%dtau(ieqn)

                                ! Need to compute row and column extends in diagonal so we can
                                ! selectively apply the mass matrix to the sub-block diagonal
                                nterms = data%mesh%domain(idom)%elems(ielem)%nterms_s
                                rstart = 1 + (ieqn-1) * nterms
                                rend   = (rstart-1) + nterms
                                cstart = rstart                 ! since it is square
                                cend   = rend                   ! since it is square

                                ! Add mass matrix divided by dt to the block diagonal
                                imat = lhs%dom(idom)%lblks(ielem,itime)%get_diagonal()
                                lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(rstart:rend,cstart:cend)  =  lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(rstart:rend,cstart:cend)  +  data%mesh%domain(idom)%elems(ielem)%mass*(ONE/dtau)

                            end do !ieqn
                        end do !itime
                    end do !ielem
                end do !idom
                endif



                !
                ! Assign rhs to b, which should allocate storage
                !
                !b = (rhs)  ! BEWARE: this causes an error. Parentheses operator not defined
                b = (-ONE)*rhs




                !
                ! We need to solve the matrix system Ax=b for the update vector x (dq)
                !

                ! Set forcing term. Converge 4 orders, or 1.e-8
                !forcing_term = resid/10000._rk
                !linear_solver%tol = max(1.e-8_rk, forcing_term)

                ! Solve system for newton step, dq
 

        do ivec = 1,size(v)
            call v(ivec)%clear()
            call z(ivec)%clear()
        end do




        !
        ! Allocate hessenberg matrix to store orthogonalization
        !
        h       = ZERO
        htmp    = ZERO
        dot_tmp = ZERO



        !
        ! Allocate vectors for solving hessenberg system
        !
        p = ZERO
        y = ZERO
        c = ZERO
        s = ZERO



        !
        ! Set initial solution x. ZERO
        !
        q = qold
        q0 = q
        x0 = dq 
        call x0%clear()
        call dq%clear()


        gmres_niter = 0


        !res = 1000000000000._rk
        res = huge(1._rk)
        do while (res > linear_solver%tol)

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
            !r0     = self%residual(A,x0,b)
            r0 =    rhs 
            r0norm = r0%norm(ChiDG_COMM)
            v(1)   = r0/r0norm
            p(1)   = r0norm



            !
            ! Inner GMRES restart loop
            !
            nvecs = 0
            do j = 1,gmres_m


                nvecs = nvecs + 1


           
                !
                ! Apply preconditioner:  z(j) = Minv * v(j)
                !
                !call timer_precon%start()
                z(j) = preconditioner%apply(lhs,v(j))
                !call timer_precon%stop()


                !
                ! Compute w = Av for the current iteration
                !
                !call timer_mv%start()
                eps = 1.0e-10_rk
                eps = 0.0_rk
                
                q = q0 + eps*z(j)
                call system%assemble( data,             &
                                      timing=timing,    &
                                      differentiate=.false.)!controller%update_lhs(lhs,niter,residual_ratio) )
                w = (rhs+b)/eps
                !w = chidg_mv(A,z(j))
                !call timer_mv%stop()






                norm_before = w%norm(ChiDG_COMM)



                !call timer_dot%start()
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
                !call timer_dot%stop()




                !call timer_norm%start()
                h(j+1,j) = w%norm(ChiDG_COMM)
                norm_after = h(j+1,j)
                !call timer_norm%stop()
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



                    !call timer_norm%start()
                    h(j+1,j) = w%norm(ChiDG_COMM)
                    !call timer_norm%stop()
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
                gmres_niter = gmres_niter + 1_ik



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
            dq = x0
            do isol = 1,nvecs
                dq = dq + y_dim(isol)*z(isol)
            end do



            !
            ! Test exit condition
            !
            if ( converged ) then
                exit
            else
                x0 = dq 
            end if



        end do   ! while







        !
        ! Report
        !
!        !err = self%error(A,x,b)
!        call self%timer%stop()
!        !call write_line('   Linear Solver Error: ',         err,                  delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<4))
!        call write_line('   Linear Solver compute time: ',  self%timer%elapsed(), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<4))
!        call write_line('   Linear Solver Iterations: ',    self%niter,           delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<4))
!
!        !call self%timer%report('Linear solver compute time: ')
!        !call timer_precon%report('Preconditioner time: ')
!        call write_line('   Preconditioner time: ',           timer_precon%elapsed(),                     delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))
!        call write_line('       Precon time per iteration: ', timer_precon%elapsed()/real(self%niter,rk), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))
!
!
!        !call timer_mv%report('MV time: ')
!        call write_line('   MV time: ',                   timer_mv%elapsed(),                     delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))
!        call write_line('       MV time per iteration: ', timer_mv%elapsed()/real(self%niter,rk), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))
!
!        !call timer_comm%report('MV comm time: ')
!        call write_line('   MV comm time: ',                   timer_comm%elapsed(),                     delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))
!        call write_line('       MV comm time per iteration: ', timer_comm%elapsed()/real(self%niter,rk), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))
!        !call timer_blas%report('MV blas time: ')
!        call write_line('   MV blas time: ',                   timer_blas%elapsed(),                     delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))
!        call write_line('       MV blas time per iteration: ', timer_blas%elapsed()/real(self%niter,rk), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))
!        call timer_comm%reset()
!        call timer_blas%reset()
!
!        !call timer_dot%report('Dot time: ')
!        !call timer_norm%report('Norm time: ')
!        call write_line('   Dot time: ',  timer_dot%elapsed(),  delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))
!        call write_line('   Norm time: ', timer_norm%elapsed(), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<5))


                !call linear_solver%solve(lhs,dq,b,preconditioner,controller,data,system)
                call self%matrix_iterations%push_back(linear_solver%niter)
                call self%matrix_time%push_back(linear_solver%timer%elapsed())




                !
                ! Line Search for appropriate step
                !
                q0 = qold
                f0 = resid
                step = 0
                if (trim(self%search) == 'Backtrack') then

                    searching = .true.
                    do while (searching)

                        !
                        ! Set line increment via backtracking.
                        !   Try: 1, 0.5, 0.25... 2^-i
                        !
                        alpha = TWO**(-real(step,rk)) 
                        call write_line("       Testing newton direction with 'alpha' = ", alpha, io_proc=GLOBAL_MASTER, silence=(verbosity<3))


                        !
                        ! Advance solution along newton direction
                        !
                        qn = q0 + alpha*dq


                        !
                        ! Clear working vector
                        !
                        !call rhs%clear()


                        !
                        ! Set working solution. Test residual at (q). Do not differentiate
                        !
                        q = qn
                        call system%assemble(data,timing=timing,differentiate=.false.)

                        !
                        ! Compute new function value
                        !
                        fn = rhs%norm(ChiDG_COMM)


                        !
                        ! Test for |R| increasing too much or NaN. 
                        !   If residual is reasonably large, still allow some growth.
                        !   If the residual is small enough, we don't want any growth.
                        ! 
                        !
                        if (ieee_is_nan(fn)) then
                            searching = .true.
                        else if ( (fn > 1.e-3_rk) .and. (fn > 2.0_rk*f0) ) then
                            searching = .true.
                        else if ( (fn < 1.e-3_rk) .and. (fn > f0) ) then
                            searching = .true.
                        else
                            searching = .false.
                        end if

                        call write_line("       Rn(Q) = ", fn, io_proc=GLOBAL_MASTER, silence=(verbosity<3))

                        step = step + 1

                    end do

                else
                    qn = q0 + dq
                end if


                !
                ! Accept new solution, clear working storage, iterate
                !
                q = qn
                call dq%clear()







                absolute_convergence = (resid > self%tol)
                relative_convergence = (fn/resid0 > self%rtol)
                !relative_convergence = ( (log10(resid0) - log10(resid)) < real(self%norders_reduction,rk) )



                ! Print iteration information
                call write_line(niter, resid, cfln(1), linear_solver%niter, update_preconditioner, update_preconditioner, delimiter='', columns=.True., column_width=30, io_proc=GLOBAL_MASTER, silence=(verbosity<2))


                call MPI_Barrier(ChiDG_COMM,ierr)
            end do ! niter


            !
            ! stop timer
            !
            call self%timer%stop()
            call self%total_time%push_back(self%timer%elapsed())
            call write_line('Nonlinear Solver elapsed time: ', self%timer%elapsed(), delimiter='', io_proc=GLOBAL_MASTER, silence=(verbosity<3))



        end associate



        call self%newton_iterations%push_back(niter)


    end subroutine solve
    !*****************************************************************************************











    !>
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/17/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute_pseudo_timestep(data,cfln)
        type(chidg_data_t),     intent(inout)   :: data
        real(rk),               intent(in)      :: cfln(:)

        integer(ik) :: idom, eqn_ID

        !
        ! Loop through elements and compute time-step function
        !
        do idom = 1,data%mesh%ndomains()

            eqn_ID = data%mesh%domain(idom)%elems(1)%eqn_ID
            call data%eqnset(eqn_ID)%compute_pseudo_timestep(idom,data%mesh,data%sdata,cfln,itime = 1)

        end do !idom

    end subroutine compute_pseudo_timestep
    !*****************************************************************************************



    !> Compute the system residual
    !!
    !! Given the system: 
    !!      Ax = b
    !!
    !!
    !! The following should be true:
    !!      b - Ax = 0
    !!
    !!
    !! So a residual is defined as:
    !!      R = b - Ax
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/17/2016
    !!
    !--------------------------------------------------------------------------------------------------------
!    function linear_residual(self,A,x,b) result(r)
!        class(jfnk_t),  intent(inout)   :: self
!        type(chidg_matrix_t),    intent(inout)   :: A
!        type(chidg_vector_t),    intent(inout)   :: x
!        type(chidg_vector_t),    intent(inout)   :: b
!
!
!        type(chidg_vector_t) :: r
!        real(rk)            :: err
!        integer(ik)         :: iparent, ielem, iblk
!
!
!        r = x
!        call r%clear()
!
!
!        !
!        ! Compute r = b - Ax
!        !
!        !r = b - A*x
!        r = b - chidg_mv(A,x)
!
!
!    end function linear_residual
!    !********************************************************************************************************









    
    subroutine destructor(self)
        type(jfnk_t),      intent(in) :: self

    end subroutine




end module type_jfnk






