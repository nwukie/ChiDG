module backward_euler_subiteration
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
    type, extends(time_scheme_t), public :: backward_euler_subiteration_t



    contains
        procedure   :: solve

        final :: destructor
    end type backward_euler_subiteration_t
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
        class(backward_euler_subiteration_t),            intent(inout)   :: self
        type(domain_t),                     intent(inout)   :: domain
        class(matrixsolver_t), optional,    intent(inout)   :: matrixsolver

        character(100)          :: filename
        integer(ik)             :: itime, nsteps, ielem, wcount, iblk, iindex, ninner, iinner, ieqn
        integer(ik)             :: rstart, rend, cstart, cend, nterms
        real(rk)                :: resid, rnorm_0, rnorm_n, dtau, amp, cfl
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

            !
            ! TIME STEP LOOP
            !
            do itime = 1,self%nsteps
                print*, "Step: ", itime


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
                    print*, "   ninner: ", ninner






                    ! Store the value of the current inner iteration solution (k) for the solution update (n+1), q_(n+1)_k
                    qold = q


                    ! Update Spatial Residual and Linearization (rhs, lin)
                    call update_space(domain)




                    ! Add mass/dt to sub-block diagonal in dR/dQ
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
                                !lin%lblks(ielem,iblk)%mat(rstart:rend,cstart:cend)  =  lin%lblks(ielem,iblk)%mat(rstart:rend,cstart:cend)  +  domain%mesh%elems(ielem)%mass/self%dt
                                !lin%lblks(ielem,iblk)%mat(rstart:rend,cstart:cend)  =  lin%lblks(ielem,iblk)%mat(rstart:rend,cstart:cend)  +  domain%mesh%elems(ielem)%mass*(ONE/dtau + ONE/self%dt)
                                lin%lblks(ielem,iblk)%mat(rstart:rend,cstart:cend)  =  lin%lblks(ielem,iblk)%mat(rstart:rend,cstart:cend)  +  domain%mesh%elems(ielem)%mass*(ONE/self%dt)
                            end if

                        end do
                    end do



                    ! Divide pseudo-time derivative by dt and multiply by mass matrix
                    dqdtau = (qold - qn)/self%dt
                    do ielem = 1,domain%mesh%nelem
                        do ieqn = 1,domain%eqnset%neqns
                            vals = matmul(domain%mesh%elems(ielem)%mass,dqdtau%lvecs(ielem)%getvar(ieqn))
                            call dqdtau%lvecs(ielem)%setvar(ieqn,vals)
                        end do
                    end do


                    ! Assign rhs to b, which should allocate storage
                    !b = (rhs)  ! BEWARE: this causes an error. Parentheses operator not defined
                    b = (-ONE)*dqdtau - rhs






                    ! We need to solve the matrix system Ax=b for the update vector x (dq)
                    call matrixsolver%solve(lin,dq,b)




                    ! Advance solution with update vector
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
                    call lin%clear()

                    


                    ! Store updated solution vector (qnew) to working solution vector (q)
                    q = qnew


                    call cpu_time(tstop)
                    telapsed = tstop - tstart
                    print*, "   Iteration time (s): ", telapsed
                    print*, "   DQ - Norm: ", resid



                    dtau = cfl * amp
                    print*, "   dtau (ps): ", dtau
                end do ! ninner

                call self%nnewton_iterations%push_back(ninner)


                if (wcount == self%nwrite) then
                    write(filename, "(I7,A4)") 1000000+itime, '.plt'
                    call write_tecio_variables(domain,trim(filename),itime+1)
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






