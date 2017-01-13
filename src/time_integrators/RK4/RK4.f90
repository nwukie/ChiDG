module RK4
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO,HALF,ONE,TWO,SIX
    use type_time_integrator,   only: time_integrator_t
    use type_chidg_data,        only: chidg_data_t
    use type_nonlinear_solver,  only: nonlinear_solver_t
    use type_linear_solver,     only: linear_solver_t
    use type_preconditioner,    only: preconditioner_t
    use type_chidgVector,

    use mod_spatial,            only: update_space
    
    use mod_tecio,              only: write_tecio_variables_unstructured

    implicit none
    private



    !> Solution advancement via the fourth order Runge-Kutta method (RK4)
    !!
    !! Given the system of partial differential equations consisting of the time-derivative
    !! of the solution vector and a spatial residual as
    !!
    !! \f$ \frac{\partial Q}{\partial t} + R(Q) = 0 \f$
    !!
    !! The time derivative is discretized as
    !!
    !! \f$ \frac{Q^{n+1} - Q^{n}}{\Delta t} - \frac{1}{6}(k_{1} + 2k_{2} + 2k_{3} + k_{4}) = 0 \f$
    !!
    !! The solution of the next time level is then computed as
    !!
    !! \f$ Q^{n+1} = Q^{n} + \Delta t \frac{1}{6}(k_{1} + 2k_{2} + 2k_{3} + k_{4}) = 0 \f$
    !! 
    !! or 
    !!
    !! \f$ Q^{n+1} = Q^{n} + \Delta Q \f$
    !!
    !! where
    !! 
    !! \f$ \Delta Q = \Delta t \frac{1}{6}(k_{1} + 2k_{2} + 2k_{3} + k_{4}) = 0 \f$
    !!
    !! with
    !!
    !! \f$ k_{1} = -inv(M)*R(Q^{n}); k_{2} = -inv(M)*R(Q^{n} + 0.5\Delta t k_{1}) \f$ and
    !! \f$ k_{3} = -inv(M)*R(Q^{n} + 0.5\Delta t k_{2}); k_{4} = -inv(M)*R(Q^{n} + \Delta t k_{3}) \f$
    !!
    !! where M refers to the mass matrix
    !!
    !!
    !! This routine computes \f$ \Delta Q \f$ and updates the solution
    !!
    !!
    !! @author Mayank Sharma
    !! @date   1/13/2017
    !!
    !!
    !----------------------------------------------------------------------------------------------------------
    type, extends(time_integrator_t), public :: RK4_t


    contains
        procedure   :: iterate

        final       :: destructor
    end type RK4_t
    !**********************************************************************************************************



    !> Solve for update dq
    !!
    !! \f$ \Delta Q = \frac{1}{6}(dQ_{1} + 2dQ_{2} + 2dQ_{3} + dQ_{4}) \f$
    !!
    !! @author Mayank Sharma
    !! @date   1/13/2017
    !!
    !----------------------------------------------------------------------------------------------------------
    subroutine iterate(self,data,nonlinear_solver,linear_solver,preconditioner)
        class(RK4_t),                           intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        class(nonlinear_solver_t),  optional,   intent(inout)   :: nonlinear_solver
        class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner

        type(chidgVector_t)     :: dq_temp(4)  ! TODO: Check if necessary
        type(chidgVector_t)     :: q_n
        character(100)          :: filename
        integer(ik)             :: itime = 1, nsteps, ielem, wcount, iblk, ieqn, idom, istage
        real(rk),allocatable    :: vals(:)


        wcount = 1
        associate( q => data%sdata%q, dq => data%sdata%dq, rhs => data%sdata%rhs, lhs => data%sdata%lhs, dt => data%sdata%dt)

            print *, 'entering time'
            do itime = 1,self%nsteps
                print *, 'Step: ', itime

                !
                ! Copy the solution vector at the present time level to a temporary vector
                !
                q_n = q

                !
                ! Update Spatial Residual and Linearization (rhs, lin)
                !
                call update_space(data)

                !
                ! Loop through the 4 increments of the fourth order Runge-Kutta method
                !
                do istage = 1,4
                    do idom = 1,data%ndomains()
                        do ielem = 1,data%mesh(idom)%nelem
                            
                            do ieqn = 1,data%eqnset(idom)%prop%nprimary_fields()
                                vals = matmul(data%mesh(idom)%elems(ielem)%invmass, rhs%dom(idom)%vecs(ielem)%getvar(ieqn,itime))
                                call rhs%dom(idom)%vecs(ielem)%setvar(ieqn,itime,vals)
                            end do 

                        end do  ! ielem
                    end do  ! idom

                    select case(istage)
                        case(1)
                                !
                                ! Compute the update for the present stage and use it to update q for the next stage
                                !
                                dq              = (-dt)*rhs
                                dq_temp(istage) = dq
                                q               = q_n + HALF*dq_1

                                !
                                ! Clear residual and linearization storage 
                                ! Compute rhs for the next stage
                                !
                                call rhs%clear()
                                call lhs%clear()
                                call update_space(data)

                        case(2)
                                !
                                ! Compute the update for the present stage and use it to update q for the next stage
                                !
                                dq              = (-dt)*rhs
                                dq_temp(istage) = dq
                                q               = q_n + HALF*dq_2

                                !
                                ! Clear residual and linearization storage
                                ! Compute rhs for the next stage
                                !
                                call rhs%clear()
                                call lhs%clear()
                                call update_space(data)

                        case(3)
                                !
                                ! Compute the update for the present stage and use it to update q for the next stage
                                !
                                dq              = (-dt)*rhs
                                dq_temp(istage) = dq
                                q               = q_n + dq_3

                                !
                                ! Clear residual and linearization storage
                                ! Compute rhs for the next stage
                                !
                                call rhs%clear()
                                call lhs%clear()
                                call update_space(data)

                        case(4)
                                !
                                ! Compute the update for the present stage
                                !
                                dq              = (-dt)*rhs
                                dq_temp(istage) = dq 

                    end select

                end do  ! istage


                !
                ! Compute update vector for the next time level
                ! Refer to the RK4 scheme description above
                !
                dq = (ONE/SIX)*(dq_temp(1) + TWO*dq_temp(2) + TWO*dq_temp(3) + dq_temp(4))


                !
                ! Advance solution with update vector
                !
                q = q_n + dq



                !
                ! Print diagnostics
                !
                call write_line('   R(Q) - Norm:    ', rhs%norm(ChiDG_COMM),delimiter='')

                
                !
                ! TODO: How does vtk fit in here?
                !
                if (wcount == self%nwrite) then
                    write(filename, '(I7,A4)') 1000000 + itime, '.plt'
                    call write_tecio_variables_unstructured(data,time(filename),itime + 1)
                    wcount = 0
                end if


                ! Clear residual and linearization storage
                call rhs%clear()
                call lhs%clear()

                wcount = wcount + 1

            end do  ! itime

        end associate


    end subroutine iterate
    !**********************************************************************************************************



    !>
    !!
    !! @author Mayank Sharma
    !! @date   1/13/2017
    !!
    !----------------------------------------------------------------------------------------------------------
    subroutine destructor(self)
        type(RK4_t),    intent(in) :: self

    end subroutine destructor
    !**********************************************************************************************************




















end module RK4
