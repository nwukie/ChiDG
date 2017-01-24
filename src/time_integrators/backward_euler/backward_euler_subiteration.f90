module backward_euler_subiteration
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, DIAG
    use type_time_scheme,       only: time_scheme_t
    use type_chidg_data,        only: chidg_data_t
    use type_nonlinear_solver,  only: nonlinear_solver_t
    use type_linear_solver,     only: linear_solver_t
    use type_preconditioner,    only: preconditioner_t
    use type_chidg_vector

    use mod_spatial,            only: update_space

    use mod_tecio,              only: write_tecio_variables
    implicit none
    private






    !>  Solution advancement via the backward-euler method
    !!
    !!  @author Nathan A. Wukie
    !!
    !----------------------------------------------------------------------------------------------------
    type, extends(time_scheme_t), public :: backward_euler_subiteration_t



    contains

        procedure   :: iterate

        final       :: destructor

    end type backward_euler_subiteration_t
    !*****************************************************************************************************








contains





    !> Solve for update 'dq'
    !!
    !!
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine iterate(self,data,nonlinear_solver,linear_solver,preconditioner)
        class(backward_euler_subiteration_t),   intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        class(nonlinear_solver_t),  optional,   intent(inout)   :: nonlinear_solver
        class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner

        character(100)          :: filename
        integer(ik)             :: itime, nsteps, ielem, wcount, iblk, iindex, ninner, iinner, ieqn
        integer(ik)             :: rstart, rend, cstart, cend, nterms, idom
        real(rk)                :: resid, rnorm_0, rnorm_n, dtau, amp, cfl
        real                    :: tstart, tstop, telapsed
        real(rk), allocatable   :: vals(:)
        type(chidg_vector_t)     :: b, qn, qold, qnew, dqdtau
      
        integer(ik)             :: ninner_iterations(self%nsteps)    ! Record number of inner iterations for each step




        tstart = 0.
        tstop = 0.
        telapsed = 0.

        wcount = 1
        ninner = 10
        associate ( q => data%sdata%q, dq => data%sdata%dq, rhs => data%sdata%rhs, lhs => data%sdata%lhs, dt => self%dt)

            call write_line('Entering time')

            !
            ! TIME STEP LOOP
            !
            do itime = 1,self%nsteps
                call write_line("Step: ", itime, delimiter='')


                ! Store qn, since q will be operated on in the inner loop
                qn = q


                !
                ! NONLINEAR CONVERGENCE INNER LOOP
                !
                resid  = ONE    ! Force inner loop entry
                ninner = 0      ! Initialize inner loop counter
                cfl    = 1._rk
                dtau   = 0.1_rk
                do while ( resid > self%tol )
                    call cpu_time(tstart)
                    ninner = ninner + 1
                    call write_line("   ninner: ", ninner)






                    ! Store the value of the current inner iteration solution (k) for the solution update (n+1), q_(n+1)_k
                    qold = q



                    !
                    ! Update Spatial Residual and Linearization (rhs, lhs)
                    !
                    call update_space(data)




                    ! Add mass/dt to sub-block diagonal in dR/dQ
                    do idom = 1,data%ndomains()
                        do ielem = 1,data%mesh(idom)%nelem
                            nterms = data%mesh(idom)%nterms_s

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
                                    lhs%dom(idom)%lblks(ielem,iblk)%mat(rstart:rend,cstart:cend)  =  lhs%dom(idom)%lblks(ielem,iblk)%mat(rstart:rend,cstart:cend)  +  data%mesh(idom)%elems(ielem)%mass*(ONE/self%dt)
                                end if

                            end do ! ieqn

                        end do ! ielem
                    end do ! idom



                    !
                    ! Divide pseudo-time derivative by dt and multiply by mass matrix
                    !
                    dqdtau = (qold - qn)/self%dt
                    do idom = 1,data%ndomains()
                        do ielem = 1,data%mesh(idom)%nelem
                            do ieqn = 1,data%eqnset(idom)%item%neqns

                                vals = matmul(data%mesh(idom)%elems(ielem)%mass,dqdtau%dom(idom)%lvecs(ielem)%getvar(ieqn))
                                call dqdtau%dom(idom)%lvecs(ielem)%setvar(ieqn,vals)

                            end do
                        end do
                    end do


                    !
                    ! Assign rhs to b, which should allocate storage
                    !
                    !b = (rhs)  ! BEWARE: this causes an error. Parentheses operator not defined
                    b = (-ONE)*dqdtau - rhs




                    !
                    ! We need to solve the matrix system Ax=b for the update vector x (dq)
                    !
                    call linear_solver%solve(lhs,dq,b,preconditioner)



                    !
                    ! Advance solution with update vector
                    !
                    qnew = qold + dq



                    ! Compute residual of nonlinear iteration
                    resid = dq%norm()

                    ! Store residual norm for first iteration
                    if (ninner == 1) then
                        rnorm_0 = dq%norm()
                    end if

                    ! Store current residual norm
                    rnorm_n = dq%norm()
                    amp = rnorm_0/rnorm_n

                    ! Clear working storage
                    call rhs%clear()
                    call dq%clear()
                    call lhs%clear()

                    

                    !
                    ! Store updated solution vector (qnew) to working solution vector (q)
                    !
                    q = qnew


                    call cpu_time(tstop)
                    telapsed = tstop - tstart
                    call write_line("   Iteration time (s): ", telapsed, delimiter='')
                    call write_line("   DQ - Norm: ", resid, delimiter='')



                end do ! ninner

                call self%newton_iterations%push_back(ninner)


                if (wcount == self%nwrite) then
                    write(filename, "(I7,A4)") 1000000+itime, '.plt'
                    call write_tecio_variables(data,trim(filename),itime+1)
                    wcount = 0
                end if
                wcount = wcount + 1


            end do  ! itime

        end associate



    end subroutine solve







    
    subroutine destructor(self)
        type(backward_euler_subiteration_t),      intent(in) :: self

    end subroutine




end module backward_euler_subiteration






