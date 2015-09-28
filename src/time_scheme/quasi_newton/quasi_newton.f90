module quasi_newton
    use mod_kinds,          only: rk,ik
    use mod_constants,      only: ZERO, ONE, TWO, DIAG
    use atype_time_scheme,  only: time_scheme_t
    use type_domain,        only: domain_t
    use atype_matrixsolver, only: matrixsolver_t
    use type_blockvector

    use mod_spatial,        only: update_space

    use mod_tecio,          only: write_tecio_variables

    use mod_entropy,        only: compute_entropy_error
    implicit none
    private



    !>  Solution advancement via the backward-euler method
    !!
    !!  @author Nathan A. Wukie
    !!
    !------------------------------------------------------------
    type, extends(time_scheme_t), public :: quasi_newton_t



    contains
        procedure   :: solve

        final :: destructor
    end type quasi_newton_t
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
        class(quasi_newton_t),              intent(inout)   :: self
        type(domain_t),                     intent(inout)   :: domain
        class(matrixsolver_t), optional,    intent(inout)   :: matrixsolver

        character(100)          :: filename
        integer(ik)             :: itime, nsteps, ielem, wcount, iblk, iindex, ninner, iinner, ieqn
        integer(ik)             :: rstart, rend, cstart, cend, nterms
        real(rk)                :: resid, rnorm_0, rnorm_n, dtau, amp, cfl, cfl0, cfln, entropy_error
        real(rk), allocatable   :: vals(:)
        type(blockvector_t)     :: b, qn, qold, qnew, dqdtau
      
        integer(ik)             :: ninner_iterations(self%nsteps)    ! Record number of inner iterations for each step





        wcount = 1
        ninner = 10
        associate ( q => domain%sdata%q, dq => domain%sdata%dq, rhs => domain%sdata%rhs, lin => domain%sdata%lin, dt => self%dt)

            print*, 'entering time'
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
            resid  = ONE    ! Force inner loop entry
            ninner = 0      ! Initialize inner loop counter
            cfl0    = 1._rk
            amp    = 0.000000000001_rk
            dtau   = cfl * amp

            cfln = cfl0



            do while ( resid > self%tol )
                ninner = ninner + 1
                print*, "   ninner: ", ninner



                !dtau = dcfln/30._rk
                dtau = dtau * 5._rk



                ! Store the value of the current inner iteration solution (k) for the solution update (n+1), q_(n+1)_k
                qold = q


                !
                ! Update Spatial Residual and Linearization (rhs, lin)
                !
                call update_space(domain)



                !
                ! Add mass/dt to sub-block diagonal in dR/dQ
                !
                do ielem = 1,domain%mesh%nelem
                    nterms = domain%mesh%nterms_s
                    do ieqn = 1,domain%eqnset%neqns
                        iblk = DIAG
                        ! Need to compute row and column extends in diagonal so we can
                        ! selectively apply the mass matrix to the sub-block diagonal
                        rstart = 1 + (ieqn-1) * nterms
                        rend   = (rstart-1) + nterms
                        cstart = rstart                 ! since it is square
                        cend   = rend                   ! since it is square
                   
                        if (allocated(lin%lblks(ielem,iblk)%mat)) then
                            ! Add mass matrix divided by dt to the block diagonal
                            lin%lblks(ielem,iblk)%mat(rstart:rend,cstart:cend)  =  lin%lblks(ielem,iblk)%mat(rstart:rend,cstart:cend)  +  domain%mesh%elems(ielem)%mass*(ONE/dtau)
                        end if

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
                call matrixsolver%solve(lin,dq,b)



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
                cfln = cfl0*(rnorm_0/rnorm_n)*10._rk

                ! Clear working storage
                call rhs%clear()
                call dq%clear()
                call lin%clear()

                


                ! Store updated solution vector (qnew) to working solution vector (q)
                q = qnew


                print*, "   DQ - Norm: ", resid
                print*, "   dtau (ps): ", dtau


                if (wcount == self%nwrite) then
                    write(filename, "(I7,A4)") 1000000+ninner, '.plt'
                    call write_tecio_variables(domain,trim(filename),ninner+1)
                    wcount = 0
                end if
                wcount = wcount + 1


            end do ! ninner

            ninner_iterations(1) = ninner   ! Record number of inner iterations

            !
            ! stop timer
            !
            call self%timer%stop()
            call self%timer%report('Solver Elapsed Time:')


            ! Write Final Solution
            write(filename, "(I7,A4)") 1000000+ninner, '.plt'
            call write_tecio_variables(domain,trim(filename),ninner+1)


        end associate



        self%ninner_iterations = ninner_iterations  ! store inner iteration count to time-scheme object



        entropy_error = compute_entropy_error(domain)
        print*, 'Entropy error: ', entropy_error

    end subroutine solve







    
    subroutine destructor(self)
        type(quasi_newton_t),      intent(in) :: self

    end subroutine




end module quasi_newton






