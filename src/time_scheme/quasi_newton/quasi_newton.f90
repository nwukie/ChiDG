module quasi_newton
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, DIAG
    use atype_time_scheme,      only: time_scheme_t
    use type_domain,            only: domain_t
    use atype_matrixsolver,     only: matrixsolver_t
    use type_preconditioner,    only: preconditioner_t
    use type_blockvector

    use mod_spatial,            only: update_space

    use mod_tecio,              only: write_tecio_variables

    use mod_entropy,            only: compute_entropy_error
    use mod_timestep,           only: compute_timestep
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
    subroutine solve(self,domain,matrixsolver,preconditioner)
        class(quasi_newton_t),                  intent(inout)   :: self
        type(domain_t),                         intent(inout)   :: domain
        class(matrixsolver_t),      optional,   intent(inout)   :: matrixsolver
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner

        character(100)          :: filename
        integer(ik)             :: itime, nsteps, ielem, wcount, iblk, iindex, ninner, iinner, ieqn
        integer(ik)             :: rstart, rend, cstart, cend, nterms
        real(rk)                :: rnorm0, rnorm, dtau, amp, cfl, cfl0, cfln, entropy_error, timing
        real(rk), allocatable   :: vals(:)
        type(blockvector_t)     :: b, qn, qold, qnew, dqdtau
      




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
            rnorm  = ONE    ! Force inner loop entry
            ninner = 0      ! Initialize inner loop counter
            !cfl0   = 2._rk




            do while ( rnorm > self%tol )
                ninner = ninner + 1
                print*, "   ninner: ", ninner



                !
                ! Store the value of the current inner iteration solution (k) for the solution update (n+1), q_(n+1)_k
                !
                qold = q


                !
                ! Update Spatial Residual and Linearization (rhs, lin)
                !
                call update_space(domain,timing)
                call self%residual_time%push_back(timing)




                !
                ! Compute and store residual norm
                !
                ! Store residual norm for first iteration
                if (ninner == 1) then
                    rnorm0 = rhs%norm()
                end if
                rnorm = rhs%norm()



                !
                ! Print and store residual
                !
                print*, "   R - Norm: ", rnorm
                call self%residual_norm%push_back(rnorm)



                !
                ! Compute new cfl for pseudo-timestep
                !
                cfln = self%cfl0*(rnorm0/rnorm)



                !
                ! Compute element-local pseudo-timestep
                !
                call compute_timestep(domain,cfln)



                !
                ! Add mass/dt to sub-block diagonal in dR/dQ
                !
                do ielem = 1,domain%mesh%nelem
                    nterms = domain%mesh%nterms_s   ! get number of solution terms
                    dtau   = domain%sdata%dt(ielem) ! get element-local timestep

                    !
                    ! Loop through equations and add mass matrix
                    !
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
                call matrixsolver%solve(lin,dq,b,preconditioner)
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
                call lin%clear()

                

                !
                ! Store updated solution vector (qnew) to working solution vector (q)
                !
                q = qnew



                if (wcount == self%nwrite) then
                    write(filename, "(I7,A4)") 1000000+ninner, '.plt'
                    call write_tecio_variables(domain,trim(filename),ninner+1)
                    wcount = 0
                end if
                wcount = wcount + 1


            end do ! ninner


            !
            ! stop timer
            !
            call self%timer%stop()
            call self%timer%report('Solver Elapsed Time:')
            call self%total_time%push_back(self%timer%elapsed())


            ! Write Final Solution
            write(filename, "(I7,A4)") 1000000+ninner, '.plt'
            call write_tecio_variables(domain,trim(filename),ninner+1)


        end associate



        call self%newton_iterations%push_back(ninner)



        entropy_error = compute_entropy_error(domain)
        print*, 'Entropy error: ', entropy_error

    end subroutine solve







    
    subroutine destructor(self)
        type(quasi_newton_t),      intent(in) :: self

    end subroutine




end module quasi_newton






