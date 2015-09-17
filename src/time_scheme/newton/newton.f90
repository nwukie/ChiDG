module newton
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO, ONE, TWO, DIAG
    use atype_time_scheme,  only: time_scheme_t
    use type_domain,        only: domain_t
    use atype_matrixsolver, only: matrixsolver_t
    use type_blockvector

    use mod_spatial,    only: update_space

    use mod_tecio,      only: write_tecio_variables
    implicit none
    private



    !>  Solution advancement via the backward-euler method
    !!
    !!  @author Nathan A. Wukie
    !!
    !------------------------------------------------------------
    type, extends(time_scheme_t), public :: newton_t



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
    subroutine solve(self,domain,matrixsolver)
        class(newton_t),              intent(inout)   :: self
        type(domain_t),                     intent(inout)   :: domain
        class(matrixsolver_t), optional,    intent(inout)   :: matrixsolver

        character(100)          :: filename
        integer(ik)             :: itime, nsteps, ielem, wcount, iblk, iindex, ninner, iinner, ieqn
        integer(ik)             :: rstart, rend, cstart, cend, nterms
        real(rk)                :: resid, rnorm_0, rnorm_n, dtau, amp, cfl, cfl0, cfln
        real                    :: tstart, tstop, telapsed
        real(rk), allocatable   :: vals(:)
        type(blockvector_t)     :: b, qn, qold, qnew, dqdtau
      
        integer(ik)             :: ninner_iterations(self%nsteps)    ! Record number of inner iterations for each step




        tstart = 0.
        tstop = 0.
        telapsed = 0.

        wcount = 1
        ninner = 10
        associate ( q => domain%sdata%q, dq => domain%sdata%dq, rhs => domain%sdata%rhs, lin => domain%sdata%lin, dt => self%dt)

            print*, 'entering time'


            ! Store qn, since q will be operated on in the inner loop
            qn = q


            !
            ! NONLINEAR CONVERGENCE INNER LOOP
            !
            resid  = ONE    ! Force inner loop entry
            ninner = 0      ! Initialize inner loop counter
            cfl0    = 1._rk
            amp    = 0.01_rk
            dtau   = cfl * amp

            cfln = cfl0

            do while ( resid > self%tol )
                call cpu_time(tstart)
                ninner = ninner + 1
                print*, "   ninner: ", ninner



                !dtau = dcfln/30._rk
                dtau = dtau * 10._rk



                ! Store the value of the current inner iteration solution (k) for the solution update (n+1), q_(n+1)_k
                qold = q


                ! Update Spatial Residual and Linearization (rhs, lin)
                call update_space(domain)

                resid = rhs%norm()


                ! Assign rhs to b, which should allocate storage
                !b = (rhs)  ! BEWARE: this causes an error. Parentheses operator not defined
                b = (-ONE)*rhs




                ! We need to solve the matrix system Ax=b for the update vector x (dq)
                call matrixsolver%solve(lin,dq,b)



                ! Advance solution with update vector
                qnew = qold + dq


!                ! Compute residual of nonlinear iteration
!                resid = dq%norm()
!
!                ! Store residual norm for first iteration
!                if (ninner == 1) then
!                    rnorm_0 = dq%norm()
!                end if
!
!                ! Store current residual norm
!                rnorm_n = dq%norm()
!                cfln = cfl0*(rnorm_0/rnorm_n)*10._rk

                ! Clear working storage
                call rhs%clear()
                call dq%clear()
                call lin%clear()

                


                ! Store updated solution vector (qnew) to working solution vector (q)
                q = qnew


                call cpu_time(tstop)
                telapsed = tstop - tstart
                print*, "   Iteration time (s): ", telapsed
!                print*, "   DQ - Norm: ", resid
                print*, "   R(Q) - Norm: ", resid
                print*, "   dtau (ps): ", dtau


                if (wcount == self%nwrite) then
                    write(filename, "(I7,A4)") 1000000+ninner, '.plt'
                    call write_tecio_variables(domain,trim(filename),ninner+1)
                    wcount = 0
                end if
                wcount = wcount + 1


            end do ! ninner

            ninner_iterations(1) = ninner   ! Record number of inner iterations





        end associate



        self%ninner_iterations = ninner_iterations  ! store inner iteration count to time-scheme object


    end subroutine solve







    
    subroutine destructor(self)
        type(newton_t),      intent(in) :: self

    end subroutine




end module newton






