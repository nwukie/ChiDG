module type_backward_euler
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, DIAG, XI_MIN, XI_MAX
    use type_time_integrator,   only: time_integrator_t
    use type_chidg_data,        only: chidg_data_t
    use type_nonlinear_solver,  only: nonlinear_solver_t
    use type_linear_solver,     only: linear_solver_t
    use type_preconditioner,    only: preconditioner_t
    use type_chidg_vector

    use mod_spatial,            only: update_space

    use mod_tecio,              only: write_tecio_variables_unstructured
    implicit none
    private



    !>  Solution advancement via the backward-euler method
    !!
    !!  @author Nathan A. Wukie
    !!
    !----------------------------------------------------------------------------------
    type, extends(time_integrator_t), public :: backward_euler_t

    contains
    
        procedure   :: iterate
        final       :: destructor

    end type backward_euler_t
    !----------------------------------------------------------------------------------










contains





    !> Solve for update 'dq'
    !!
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine iterate(self,data,nonlinear_solver,linear_solver,preconditioner)
        class(backward_euler_t),                intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        class(nonlinear_solver_t),  optional,   intent(inout)   :: nonlinear_solver
        class(linear_solver_t),     optional,   intent(inout)   :: linear_solver 
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner

        character(100)          :: filename
        integer(ik)             :: itime, nsteps, ielem, wcount, iblk, iindex, &
                                   ninner, iinner, ieqn, idom
        integer(ik)             :: rstart, rend, cstart, cend, nterms
        real(rk)                :: resid, rnorm_0, rnorm_n
        real                    :: tstart, tstop, telapsed
        real(rk), allocatable   :: vals(:)
        type(chidg_vector_t)     :: b, qn, qold, qnew, dtau
      


        tstart = 0.
        tstop = 0.
        telapsed = 0.

        print*, 'hi - 1'

        wcount = 1
        ninner = 10
        associate ( q => data%sdata%q, dq => data%sdata%dq, rhs => data%sdata%rhs, lhs => data%sdata%lhs, dt => self%dt)

        !
        ! TIME STEP LOOP
        !
        print*, 'hi - 2'
        do itime = 1,self%nsteps
            print*, "Step: ", itime

        print*, 'hi - 3'
            if (itime == 1) then
                write(filename, "(I7,A4)") 1000000, '.plt'
                call write_tecio_variables_unstructured(data,trim(filename),1)
            end if

        print*, 'hi - 4'
            !
            ! Store the value of the current inner iteration solution (k) for the solution update (n+1), q_(n+1)_k
            !
            qold = q


        print*, 'hi - 5'
            !
            ! Write lhs and rhs to file
            !
            open(unit=16, file="Q0.txt")

        print*, 'hi - 6'
            write(16,*) "double Q0(100) = { "
            do idom = 1,data%ndomains()
                do ielem = 1,data%mesh(idom)%nelem
                    write(16,*) q%dom(idom)%vecs(ielem)%vec(1)
                end do
            end do

            write(16,*) " } "
            close(16)




        print*, 'hi - 7'








            !
            ! Update Spatial Residual and Linearization (rhs, lin)
            !
            call update_space(data)


        print*, 'hi - 8'

            !
            ! Compute residual of nonlinear iteration
            !
            resid = rhs%norm()

            print*, "||R||: ", resid


        print*, 'hi - 9'

            !
            ! Add mass/dt to sub-block diagonal in dR/dQ
            !
            do idom = 1,data%ndomains()
                do ielem = 1,data%mesh(idom)%nelem
                    nterms = data%mesh(idom)%nterms_s
                    do ieqn = 1,data%eqnset(idom)%prop%nprimary_fields()
                        iblk = DIAG
                        ! Need to compute row and column extends in diagonal so we can
                        ! selectively apply the mass matrix to the sub-block diagonal
                        rstart = 1 + (ieqn-1) * nterms
                        rend   = (rstart-1) + nterms
                        cstart = rstart                 ! since it is square
                        cend   = rend                   ! since it is square
                   
                        if (allocated(lhs%dom(idom)%lblks(ielem,iblk)%mat)) then
                            ! Add mass matrix divided by dt to the block diagonal
                            lhs%dom(idom)%lblks(ielem,iblk)%mat(rstart:rend,cstart:cend)  =  lhs%dom(idom)%lblks(ielem,iblk)%mat(rstart:rend,cstart:cend)  +   data%mesh(idom)%elems(ielem)%mass/self%dt
                        end if

                    end do
                end do
            end do

        print*, 'hi - 10'
            !
            ! Assign rhs to b, which should allocate storage
            !
            b = (-ONE)*rhs

            !
            ! Write lhs and rhs to file
            !
            open(unit=11, file="LHS_diagonal.txt")
            open(unit=12, file="LHS_lower.txt")
            open(unit=13, file="LHS_upper.txt")
            open(unit=14, file="RHS.txt")

            write(11,*) "double LHS_diagonal(100) = { "
            write(12,*) "double LHS_lower(99) = { "
            write(13,*) "double LHS_upper(99) = { "
            write(14,*) "double RHS(100) = { "
            do idom = 1,data%ndomains()
                do ielem = 1,data%mesh(idom)%nelem
                   
                    if (allocated(lhs%dom(idom)%lblks(ielem,DIAG)%mat)) then
                        write(11,*) lhs%dom(idom)%lblks(ielem,DIAG)%mat(1,1)
                    end if

                    if (allocated(lhs%dom(idom)%lblks(ielem,XI_MIN)%mat)) then
                        write(12,*) lhs%dom(idom)%lblks(ielem,XI_MIN)%mat(1,1)
                    end if

                    if (allocated(lhs%dom(idom)%lblks(ielem,XI_MAX)%mat)) then
                        write(13,*) lhs%dom(idom)%lblks(ielem,XI_MAX)%mat(1,1)
                    end if

                    write(14,*) b%dom(idom)%vecs(ielem)%vec(1)

                end do
            end do

            write(11,*) " } "
            write(12,*) " } "
            write(13,*) " } "
            write(14,*) " } "
            close(11)
            close(12)
            close(13)
            close(14)



        print*, 'hi - 11'



            !
            ! Solve the matrix system Ax=b for the update vector x (dq)
            !
            call linear_solver%solve(lhs,dq,b,preconditioner)





            !
            ! Write lhs and rhs to file
            !
            open(unit=15, file="DQ.txt")

            write(15,*) "double DQ(100) = { "
            do idom = 1,data%ndomains()
                do ielem = 1,data%mesh(idom)%nelem
                   
                    write(15,*) dq%dom(idom)%vecs(ielem)%vec(1)

                end do
            end do

            write(15,*) " } "
            close(15)



        print*, 'hi - 12'











            !
            ! Advance solution with update vector
            !
            qnew = qold + dq


            !
            ! Compute residual of nonlinear iteration
            !
            resid = dq%norm()

        print*, 'hi - 13'

            ! Clear working storage
            call rhs%clear()
            call dq%clear()
            call lhs%clear()

            


            ! Store updated solution vector (qnew) to working solution vector (q)
            q = qnew


        print*, 'hi - 14'
            call cpu_time(tstop)
            telapsed = tstop - tstart
            print*, "   Iteration time (s): ", telapsed
            print*, "   DQ - Norm: ", resid



        print*, 'hi - 15'

            if (wcount == self%nwrite) then
                write(filename, "(I7,A4)") 1000000+itime, '.plt'
                call write_tecio_variables_unstructured(data,trim(filename),itime+1)
                wcount = 0
            end if
            wcount = wcount + 1


        end do  ! itime

        end associate





    end subroutine iterate
    !*****************************************************************************************







    
    subroutine destructor(self)
        type(backward_euler_t),      intent(in) :: self

    end subroutine




end module type_backward_euler
