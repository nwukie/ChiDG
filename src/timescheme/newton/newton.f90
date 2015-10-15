module newton
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, DIAG
    use type_timescheme,        only: timescheme_t
    use type_chidgData,         only: chidgData_t
    use atype_matrixsolver,     only: matrixsolver_t
    use type_preconditioner,    only: preconditioner_t
    use type_blockvector

    use mod_spatial,    only: update_space

    use mod_tecio,      only: write_tecio_variables

    !use mod_entropy,    only: compute_entropy_error
    implicit none
    private



    !>  Solution advancement via the newton's method
    !!
    !!
    !!
    !!
    !!
    !!  @author Nathan A. Wukie
    !!
    !------------------------------------------------------------
    type, extends(timescheme_t), public :: newton_t



    contains
        procedure   :: solve

        final :: destructor
    end type newton_t
    !-----------------------------------------------------------










contains





    !> Solve for update 'dq'
    !!
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    !subroutine solve(self,domains,matrixsolver,preconditioner)
    subroutine solve(self,data,matrixsolver,preconditioner)
        class(newton_t),                    intent(inout)   :: self
        !type(domain_t),                     intent(inout)   :: domains(:)
        type(chidgData_t),                  intent(inout)   :: data
        class(matrixsolver_t),   optional,  intent(inout)   :: matrixsolver
        class(preconditioner_t), optional,  intent(inout)   :: preconditioner

        character(100)          :: filename
        integer(ik)             :: itime, nsteps, ielem, wcount, iblk, iindex, niter, iinner, ieqn
        integer(ik)             :: rstart, rend, cstart, cend, nterms
        real(rk)                :: resid, timing
        real(rk), allocatable   :: vals(:)
        type(blockvector_t)     :: b, qn, qold, qnew
      

        real(rk)                :: entropy_error



!
!        wcount = 1
!        associate ( q => domain%sdata%q, dq => domain%sdata%dq, rhs => domain%sdata%rhs, lin => domain%sdata%lin, dt => self%dt)
!
!            print*, 'Entering time'
!            !
!            ! start timer
!            !
!            call self%timer%start()
!
!
!            ! Store qn, since q will be operated on in the inner loop
!            qn = q
!
!
!            !
!            ! NONLINEAR CONVERGENCE LOOP
!            !
!            resid  = ONE    ! Force inner loop entry
!            niter = 0       ! Initialize inner loop counter
!
!
!            do while ( resid > self%tol )
!                niter = niter + 1
!                print*, "   niter: ", niter
!
!
!                !
!                ! Store the value of the current inner iteration solution (k) for the solution update (n+1), q_(n+1)_k
!                !
!                qold = q
!
!
!                !
!                ! Update Spatial Residual and Linearization (rhs, lin)
!                !
!                call update_space(domain,timing)
!                call self%residual_time%push_back(timing)   ! non-essential record-keeping
!
!                resid = rhs%norm()
!
!
!                !
!                ! Print diagnostics
!                !
!                print*, "   R(Q) - Norm: ", resid
!                call self%residual_norm%push_back(resid)
!
!
!
!
!                !
!                ! Assign rhs to b, which should allocate storage
!                !
!                !b = (rhs)  ! BEWARE: this causes an error. Parentheses operator not defined
!                b = (-ONE)*rhs
!
!
!
!
!                !
!                ! We need to solve the matrix system Ax=b for the update vector x (dq)
!                !
!                !call matrixsolver%solve(lin,dq,b,preconditioner)
!
!                call matrixsolver%solve(domains,preconditioner)
!
!                call self%matrix_iterations%push_back(matrixsolver%niter)       ! non-essential record-keeping
!                call self%matrix_time%push_back(matrixsolver%timer%elapsed())   ! non-essential record-keeping
!
!
!
!                !
!                ! Advance solution with update vector
!                !
!                qnew = qold + dq
!
!
!
!                !
!                ! Clear working storage
!                !
!                call rhs%clear()
!                call dq%clear()
!                call lin%clear()
!
!                
!
!
!                !
!                ! Store updated solution vector (qnew) to working solution vector (q)
!                !
!                q = qnew
!
!
!
!
!
!    
!                !
!                ! Write incremental solution
!                !
!                if (wcount == self%nwrite) then
!                    write(filename, "(I7,A4)") 1000000+niter, '.plt'
!                    call write_tecio_variables(domain,trim(filename),niter+1)
!                    wcount = 0
!                end if
!                wcount = wcount + 1
!
!
!            end do ! niter
!
!
!            !
!            ! stop timer. Record timings
!            !
!            call self%timer%stop()
!            call self%timer%report('Solver Elapsed Time:')
!            call self%total_time%push_back(self%timer%elapsed())
!
!
!
!
!            !
!            ! Write final solution
!            !
!            !write(filename, "(I7,A4)") 1000000+niter, '.plt'
!            !call write_tecio_variables(domain,trim(filename),niter+1)
!
!
!
!        end associate
!
!
!
!        !
!        ! Store newton iteration count
!        !
!        call self%newton_iterations%push_back(niter)
!
!
!        entropy_error = compute_entropy_error(domain) 
!        print*, 'Entropy error: ', entropy_error
    end subroutine solve







    
    subroutine destructor(self)
        type(newton_t),      intent(in) :: self

    end subroutine




end module newton






