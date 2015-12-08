module newton
    use messenger,              only: write_line
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, DIAG
    use type_timescheme,        only: timescheme_t
    use type_chidg_data,        only: chidg_data_t
    use atype_matrixsolver,     only: matrixsolver_t
    use type_preconditioner,    only: preconditioner_t
    use type_chidgVector

    use mod_spatial,    only: update_space

    use mod_tecio,      only: write_tecio_variables

    use mod_entropy,    only: compute_entropy_error
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
    subroutine solve(self,data,matrixsolver,preconditioner)
        class(newton_t),                    intent(inout)   :: self
        type(chidg_data_t),                 intent(inout)   :: data
        class(matrixsolver_t),   optional,  intent(inout)   :: matrixsolver
        class(preconditioner_t), optional,  intent(inout)   :: preconditioner

        character(100)          :: filename
        integer(ik)             :: itime, nsteps, ielem, wcount, iblk, iindex, niter, iinner, ieqn
        integer(ik)             :: rstart, rend, cstart, cend, nterms
        real(rk)                :: resid, timing
        real(rk), allocatable   :: vals(:)
        type(chidgVector_t)     :: b, qn, qold, qnew
      

        real(rk)                :: entropy_error




        wcount = 1
        associate ( q => data%sdata%q, dq => data%sdata%dq, rhs => data%sdata%rhs, lhs => data%sdata%lhs)

            call write_line('Entering time')

            !
            ! start timer
            !
            call self%timer%start()


            ! Store qn, since q will be operated on in the inner loop
            qn = q


            !
            ! NONLINEAR CONVERGENCE LOOP
            !
            resid  = ONE    ! Force inner loop entry
            niter = 0       ! Initialize inner loop counter


            do while ( resid > self%tol )
                niter = niter + 1

                call write_line("   niter: ", niter, delimiter='')


                !
                ! Store the value of the current inner iteration solution (k) for the solution update (n+1), q_(n+1)_k
                !
                qold = q


                !
                ! Update Spatial Residual and Linearization (rhs, lin)
                !
                call update_space(data,timing)







                call self%residual_time%push_back(timing)   ! non-essential record-keeping

                resid = rhs%norm()


                !
                ! Print diagnostics
                !
                call write_line("   R(Q) - Norm: ", resid, delimiter='')
                call self%residual_norm%push_back(resid)




                !
                ! Assign rhs to b, which should allocate storage
                !
                !b = (rhs)  ! BEWARE: this causes an error. Parentheses operator not defined
                b = (-ONE)*rhs




                !
                ! We need to solve the matrix system Ax=b for the update vector x (dq)
                !
                call matrixsolver%solve(lhs,dq,b,preconditioner)

                call self%matrix_iterations%push_back(matrixsolver%niter)
                call self%matrix_time%push_back(matrixsolver%timer%elapsed())



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



                !
                ! Write incremental solution
                !
                if (wcount == self%nwrite) then
                    write(filename, "(I7,A4)") 1000000+niter, '.plt'
                    call write_tecio_variables(data,trim(filename),niter+1)
                    wcount = 0
                end if
                wcount = wcount + 1


            end do ! niter


            !
            ! stop timer. Record timings
            !
            call self%timer%stop()
            call self%timer%report('Solver Elapsed Time:')
            call self%total_time%push_back(self%timer%elapsed())



        end associate



        !
        ! Store newton iteration count
        !
        call self%newton_iterations%push_back(niter)


        entropy_error = compute_entropy_error(data)
        call write_line('Entropy error: ', entropy_error, delimiter='')

    end subroutine solve







    
    subroutine destructor(self)
        type(newton_t),      intent(in) :: self

    end subroutine




end module newton






