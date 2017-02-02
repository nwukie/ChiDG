module type_quasi_newton
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO, ONE, TWO, DIAG
    use mod_spatial,            only: update_space
    use mod_hdfio,              only: write_solution_hdf
    use mod_tecio,              only: write_tecio_variables_unstructured
    use mod_chidg_mpi,          only: ChiDG_COMM, GLOBAL_MASTER, IRANK, NRANK
    use mpi_f08,                only: MPI_Barrier

    use type_chidg_data,        only: chidg_data_t
    use type_nonlinear_solver,  only: nonlinear_solver_t
    use type_linear_solver,     only: linear_solver_t
    use type_preconditioner,    only: preconditioner_t
    use type_chidg_vector
    implicit none
    private



    !>  Solution advancement via the backward-euler method
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !------------------------------------------------------------------------------------------
    type, extends(nonlinear_solver_t), public :: quasi_newton_t


    contains
    
        procedure   :: solve
        final       :: destructor

    end type quasi_newton_t
    !******************************************************************************************










contains





    !>  Solve for update 'dq'
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine solve(self,data,linear_solver,preconditioner)
        class(quasi_newton_t),                  intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner

        character(100)          :: filename
        integer(ik)             :: itime, nsteps, ielem, wcount, iblk, iindex,  &
                                   niter, ieqn, idom, ierr,                     &
                                   rstart, rend, cstart, cend, nterms, imat, iwrite

        real(rk)                :: dtau, amp, cfl, timing, resid, resid_new
        real(rk), allocatable   :: vals(:), cfln(:), rnorm0(:), rnorm(:)
        type(chidg_vector_t)     :: b, qn, qold, qnew, dqdtau
      

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
            resid = ONE    ! Force inner loop entry
            niter = 0      ! Initialize inner loop counter


            do while ( resid > self%tol )
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
                resid = rhs%norm(ChiDG_COMM)


                !
                ! Print diagnostics, check tolerance.
                !
                call write_line("   R(Q) - Norm: ", resid, delimiter='', columns=.True., column_width=20, io_proc=GLOBAL_MASTER)
                if ( resid < self%tol ) exit
                call self%residual_time%push_back(timing)
                call self%residual_norm%push_back(resid)


                !
                ! Compute and store residual norm for each field
                !
                if (niter == 1) then
                    rnorm0 = rhs%norm_fields(ChiDG_COMM)
                end if


                rnorm = rhs%norm_fields(ChiDG_COMM)


                ! Update initial residual norm if it has increased
                where (rnorm > rnorm0)
                    rnorm0 = rnorm0 + 0.3_rk*(rnorm - rnorm0)
                end where



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


                !
                ! Compute element-local pseudo-timestep
                !
                call compute_pseudo_timestep(data,cfln)


                !
                ! Add mass/dt to sub-block diagonal in dR/dQ
                !
                do idom = 1,data%ndomains()
                    do ielem = 1,data%mesh(idom)%nelem
                        do itime = 1,data%mesh(idom)%ntime
                            do ieqn = 1,data%eqnset(idom)%prop%nprimary_fields()

                                !nterms = data%mesh(idom)%nterms_s   ! get number of solution terms
                                !dtau   = data%sdata%dt(idom,ielem)  ! get element-local timestep
                                dtau = data%mesh(idom)%elems(ielem)%dtau(ieqn)

                                ! Need to compute row and column extends in diagonal so we can
                                ! selectively apply the mass matrix to the sub-block diagonal
                                nterms = data%mesh(idom)%elems(ielem)%nterms_s
                                rstart = 1 + (ieqn-1) * nterms
                                rend   = (rstart-1) + nterms
                                cstart = rstart                 ! since it is square
                                cend   = rend                   ! since it is square

                                ! Add mass matrix divided by dt to the block diagonal
                                imat = lhs%dom(idom)%lblks(ielem,itime)%get_diagonal()
                                lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(rstart:rend,cstart:cend)  =  lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(rstart:rend,cstart:cend)  +  data%mesh(idom)%elems(ielem)%mass*(ONE/dtau)

                            end do !ieqn
                        end do !itime
                    end do !ielem
                end do !idom



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



                !
                ! Write solution if the count is right
                !
                !if (wcount == self%nwrite) then
                    if (data%eqnset(1)%get_name() == 'Navier Stokes AV') then
                        call write_solution_hdf(data,'aachen_stator_cascade.h5')
                !        write(filename,'(I2)') niter
                !        call write_tecio_variables_unstructured(data,trim(filename)//'.dat',niter)
                !        wcount = 0
                    end if
                !end if
                wcount = wcount + 1

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

        integer(ik) :: idom

        !
        ! Loop through elements and compute time-step function
        !
        do idom = 1,data%ndomains()

            call data%eqnset(idom)%compute_pseudo_timestep(idom,data%mesh,data%sdata,cfln,itime = 1)

        end do !idom

    end subroutine compute_pseudo_timestep
    !*****************************************************************************************










    
    subroutine destructor(self)
        type(quasi_newton_t),      intent(in) :: self

    end subroutine




end module type_quasi_newton






