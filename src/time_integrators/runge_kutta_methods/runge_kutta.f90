module runga_kutta
    use mod_kinds,                  only: rk,ik
    use mod_constants,              only: ZERO,HALF,ONE,TWO,SIX
    use type_time_integrator,       only: time_integrator_t
    use type_chidg_data,            only: chidg_data_t
    use type_nonlinear_solver,      only: nonlinear_solver_t
    use type_linear_solver,         only: linear_solver_t
    use type_preconditioner,        only: preconditioner_t
    use type_chidgVector,           

    use mod_spatial,                only: update_space
    
    use mod_tecio,                  only: write_tecio_variables_unstructured

    use mod_explicit_RK_methods,    
    
    implicit none
    private


    !>  Solution advancement via a general explicit Runge-Kutta method
    !!
    !!  Given the system of partial differential equations consisting of the time-derivative
    !!  of the solution vector and a spatial residual as
    !!
    !!  \f$ \frac{\partial Q}{\partial t} + R(Q) = 0 \f$
    !!
    !!  The solution of the next time level is then computed as
    !!
    !!  \f$ Q^{n+1} = Q^{n} + \sum_{i = 1}^{s} b_{i} \Delta Q_{i} \f$
    !!
    !!  or
    !!
    !!  \f$ Q^{n+1} = Q^{n} + \Delta Q \f$
    !!
    !!  where
    !!
    !!  \f$ \Delta Q = \sum_{i = 1}^{s} b_{i} \Delta Q_{i} \f$
    !!
    !!  with
    !!
    !!  \f$ \Delta Q_{i} = -\Delta t m^{-1} R(Q^{n} + \sum_{j = 1}^{i-1} a_{ij} \Delta Q_{j}) \f$
    !!
    !!  where \f$ M \f$ refers to the mass matrix, \f$ s \f$ refers to the number of the stages in the 
    !!  method and \f$ a_{ij}'s & b_{j}'s \f$ refer to the coefficents defined by the Butcher tableau
    !!
    !!
    !!  This routine computes \f$ \Delta Q \f$ and updates the solution
    !!
    !!
    !!  @author Mayank Sharma
    !!  @date   1/18/2017
    !!
    !!
    !---------------------------------------------------------------------------------------------------------------------
    type. extends(time_integrator_t), public :: runge_kutta_t


    contains
        procedure   :: iterate

        final       :: destructor
    end type runge_kutta_t
    !*********************************************************************************************************************



    !>  Solve for update dq
    !!
    !!  \f$ \Delta Q_{i} = -\Delta t m^{-1} R(Q^{n} + \sum_{j = 1}^{i-1} a_{ij} \Delta Q_{j}) \f$
    !!
    !!  @author Mayank Sharma
    !!  @date   18/1/2017
    !!
    !!
    !---------------------------------------------------------------------------------------------------------------------
    subroutine iterate(self,data,nonlinear_solver,linear_solver,preconditioner)
        class(runge_kutta_t),                   intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        class(nonlinear_solver_t),  optional,   intent(inout)   :: nonlinear_solver
        class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner

        integer(ik),                        :: nstage   ! Number of stages in the RK method
        real(rk),           allocatable     :: a(:,:)   ! a & b are coefficient arrays
        real(rk),           allocatable     :: b(:)
        type(chidgVector_t),allocatable     :: delq(:)
        type(chidgVector_t)                 :: q_n, dq_new_rhs ! dq_new_rhs is used to compute q for 
                                                               ! update_space (to compute residual for next stage)                                                                        
        character(100)                      :: filename
        integer(ik)                         :: itime = 1, nsteps, ielem, wcount, iblk, ieqn, idom, istage, j
        real(rk),           allocatable     :: vals(:)


        wcount = 1
        associate( q => data%sdata%q, dq => data%sdata%dq, rhs => data%sdata%rhs, lhs => data%sdata%lhs, dt => data%sdata%dt )
        
        !
        ! Get nstage, a and b for a particular scheme
        ! All methods defined in mod_explicit_RK_methods 
        ! TODO: Method for selecting a particular method
        !
        call runge_kutta_4th_order(nstage,a,b)

            print *, 'entering time'
            do itime = 1,self%nsteps
                print *, 'Step: ', itime

                !
                ! Copy the solution vector at the present time level to a temporary vector
                !
                q_n = q

                if (allocated(delq)) deallocate(delq)
                allocate(delq(nstage))

                !
                ! Loop through the stages of the present RK method
                !
                do istage = 1,nstage
                    
                    !
                    ! Update Spatial Residual and Linearization (rhs, lin)
                    !
                    call update_space(data)

                    do idom = 1,data%ndomains()
                        do ielem = 1,data%mesh(idom)%nelem

                            do ieqn = 1,data%eqnset(idom)%prop%nprimary_fields()
                                vals = matmul(data%mesh(idom)%elems(ielem)%invmass, rhs%dom(idom)%vecs(ielem)%getvar(ieqn,itime))
                                call rhs%dom(idom)%vecs(ielem)%setvar(ieqn,itime,vals)
                            end do

                        end do  ! ielem
                    end do  ! idom

                    !
                    ! Calculate stagewise update and its contribution to the final update
                    !
                    dq = (-dt)*rhs
                    delq(istage) = b(istage)*dq
                    
                    !
                    ! Initialize dq_new_rhs for the present stage
                    ! Computed using an accumulated sum
                    !
                    dq_new_rhs = ZERO

                    !
                    ! Compute the solution vector to pass to update_space for the next stage
                    ! Refer to scheme description above
                    !
                    if (istage < nstage) then
                        do j = 1,istage
                            dq_new_rhs = dq_new_rhs + (a(istage + 1,j)*delq(j)/b(j))
                        end do

                        q = q_n + dq_new_rhs

                    end if

                    !
                    ! Clear Residual and Linearization storage
                    !
                    call lhs%clear()
                    call rhs%clear()

                end do  ! istage


                !
                ! Compute update vector for the next time level and advance solution with update vector
                ! 
                dq = sum(delq)
                q = q_n + dq



                !
                ! Print diagnostics 
                !
                call write_line('   R(Q) - Norm:    ', rhs%norm(ChiDG_COMM),delimiter = '')


                if (wcount == self%nwrite) then
                    write(filename,'(I7,A4)') 1000000 + itime, '.plt'
                    call write_tecio_variables_unstructured(data,time(filename),itime + 1)
                    wcount = 0
                end if


                wcount = wcount + 1

            end do  ! itime

        end associate


    end subroutine iterate
    !*********************************************************************************************************************



    !>
    !!
    !!  @author Mayank Sharma
    !!  @date   1/18/2017
    !!
    !!
    !---------------------------------------------------------------------------------------------------------------------
    subroutine destructor(self)
        type(runge_kutta_t),    intent(in)  :: self

    end subroutine destructor
    !*********************************************************************************************************************




















end module runge_kutta
