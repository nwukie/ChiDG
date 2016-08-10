module quasi_newton
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, DIAG
    use mod_spatial,            only: update_space
    use mod_tecio,              only: write_tecio_variables
    use mod_chidg_mpi,          only: ChiDG_COMM, GLOBAL_MASTER

    use type_chidg_data,        only: chidg_data_t
    use type_nonlinear_solver,  only: nonlinear_solver_t
    use type_linear_solver,     only: linear_solver_t
    use type_preconditioner,    only: preconditioner_t
    use type_chidgVector



    use mod_entropy,            only: compute_entropy_error
    use mod_timestep,           only: compute_timestep
    implicit none
    private



    !>  Solution advancement via the backward-euler method
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !---------------------------------------------------------------------------------------------
    type, extends(nonlinear_solver_t), public :: quasi_newton_t


    contains
    
        procedure   :: solve
        final       :: destructor

    end type quasi_newton_t
    !**********************************************************************************************










contains





    !> Solve for update 'dq'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine solve(self,data,linear_solver,preconditioner)
        class(quasi_newton_t),                  intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner

        character(100)          :: filename
        integer(ik)             :: itime, nsteps, ielem, wcount, iblk, iindex, niter, ieqn, idom, ierr
        integer(ik)             :: rstart, rend, cstart, cend, nterms
        real(rk)                :: dtau, amp, cfl, cfln, entropy_error, timing
        real(rk)                :: rnorm0, rnorm, resid, resid_new
        real(rk), allocatable   :: vals(:)
        type(chidgVector_t)     :: b, qn, qold, qnew, dqdtau
      




        wcount = 1
        associate ( q => data%sdata%q, dq => data%sdata%dq, rhs => data%sdata%rhs, lhs => data%sdata%lhs)

            call write_line('Entering time', io_proc=GLOBAL_MASTER)
            call write_line(' ', io_proc=GLOBAL_MASTER)

            !
            ! start timer
            !
            call self%timer%reset()
            call self%timer%start()


            ! Store qn, since q will be operated on in the inner loop
            qn = q


            !
            ! NONLINEAR CONVERGENCE INNER LOOP
            !
            rnorm = ONE    ! Force inner loop entry
            niter = 0      ! Initialize inner loop counter

            dtau = 1._rk



            do while ( rnorm > self%tol )
                niter = niter + 1
                call write_line("   niter: ", niter, delimiter='', columns=.True., column_width=20, io_proc=GLOBAL_MASTER)



                !
                ! Store the value of the current inner iteration solution (k) for the solution update (n+1), q_(n+1)_k
                !
                qold = q


                !
                ! Update Spatial Residual and Linearization (rhs, lin)
                !
                call update_space(data,timing)
                call self%residual_time%push_back(timing)

                resid = rhs%norm(ChiDG_COMM)



                !
                ! Print diagnostics
                !
                call write_line("   R(Q) - Norm: ", resid, delimiter='', columns=.True., column_width=20, io_proc=GLOBAL_MASTER)
                call self%residual_norm%push_back(resid)


                !
                ! Compute and store residual norm
                !
                ! Store residual norm for first iteration
                if (niter == 1) then
                    rnorm0 = rhs%norm(ChiDG_COMM)
                end if
                rnorm = rhs%norm(ChiDG_COMM)




                !
                ! Compute new cfl for pseudo-timestep
                !
                cfln = self%cfl0*(rnorm0/rnorm)




                !
                ! Compute element-local pseudo-timestep
                !
                !data%sdata%dt = cfln
                call compute_timestep(data,cfln)




                !
                ! Add mass/dt to sub-block diagonal in dR/dQ
                !
                do idom = 1,data%ndomains()
                    do ielem = 1,data%mesh(idom)%nelem
                        nterms = data%mesh(idom)%nterms_s   ! get number of solution terms


                        dtau   = data%sdata%dt(idom,ielem)  ! get element-local timestep

                        !
                        ! Loop through equations and add mass matrix
                        !
                        do ieqn = 1,data%eqnset(idom)%item%neqns
                            iblk = DIAG
                            ! Need to compute row and column extends in diagonal so we can
                            ! selectively apply the mass matrix to the sub-block diagonal
                            rstart = 1 + (ieqn-1) * nterms
                            rend   = (rstart-1) + nterms
                            cstart = rstart                 ! since it is square
                            cend   = rend                   ! since it is square
                       
                            if (allocated(lhs%dom(idom)%lblks(ielem,iblk)%mat)) then
                                ! Add mass matrix divided by dt to the block diagonal
                                lhs%dom(idom)%lblks(ielem,iblk)%mat(rstart:rend,cstart:cend)  =  lhs%dom(idom)%lblks(ielem,iblk)%mat(rstart:rend,cstart:cend)  +  data%mesh(idom)%elems(ielem)%mass*(ONE/dtau)
                            end if

                        end do
                    end do
                end do


                !
                ! Assign rhs to b, which should allocate storage
                !
                !b = (rhs)  ! BEWARE: this causes an error. Parentheses operator not defined
                b = (-ONE)*rhs




                !
                ! We need to solve the matrix system Ax=b for the update vector x (dq)
                !
                call linear_solver%solve(lhs,dq,b,preconditioner)
                call self%matrix_iterations%push_back(linear_solver%niter)
                call self%matrix_time%push_back(linear_solver%timer%elapsed())



                !
                ! Advance solution with update vector
                !
                qnew = qold + dq



                !
                ! Clear working storage
                !
                call rhs%clear()
                call dq%clear()
                call lhs%clear()

                

                !
                ! Store updated solution vector (qnew) to working solution vector (q)
                !
                q = qnew



                !if (wcount == self%nwrite) then
                !    write(filename, "(I7,A4)") 1000000+niter, '.plt'
                !    call write_tecio_variables(data,trim(filename),niter+1)
                !    wcount = 0
                !end if
                !wcount = wcount + 1

                call MPI_Barrier(ChiDG_COMM,ierr)

            end do ! niter


            !
            ! stop timer
            !
            call self%timer%stop()
            call self%timer%report('Solver Elapsed Time: ')
            call self%total_time%push_back(self%timer%elapsed())



        end associate



        call self%newton_iterations%push_back(niter)



        !entropy_error = compute_entropy_error(data)
        !call write_line('Entropy error: ', entropy_error, delimiter='', columns=.True., column_width=20)

    end subroutine solve
    !************************************************************************************************************







    
    subroutine destructor(self)
        type(quasi_newton_t),      intent(in) :: self

    end subroutine




end module quasi_newton






