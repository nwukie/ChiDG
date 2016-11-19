module forward_euler
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: ZERO
    use type_time_integrator,   only: time_integrator_t
    use type_chidg_data,        only: chidg_data_t
    use type_nonlinear_solver,  only: nonlinear_solver_t
    use type_linear_solver,     only: linear_solver_t
    use type_preconditioner,    only: preconditioner_t
    use type_chidgVector

    use mod_spatial,            only: update_space

    use mod_tecio,              only: write_tecio_variables_unstructured
    implicit none
    private



    !>  Solution advancement via the forward-euler method
    !!
    !! Given the system of partial differential equations consisting of the time-derivative
    !! of the solution vector and a spatial residual as
    !!
    !! \f$ \frac{\partial Q}{\partial t} + R(Q) = 0 \f$
    !!
    !! the time derivative is discretized by a one-sided finite-difference approximation as
    !!
    !! \f$ \frac{Q^{n+1} - Q^{n}}{\Delta t} + R(Q^n) = 0 \f$
    !!
    !! The solution at the next time level is then computed as
    !!
    !! \f$ Q^{n+1} = Q^{n} - \Delta t R(Q^n) \f$
    !!
    !! or
    !!
    !! \f$ Q^{n+1} = Q^{n} + \Delta Q \f$
    !!
    !! where \f$ \Delta Q \f$ is defined as
    !!
    !! \f$ \Delta Q = -\Delta t R(Q) \f$
    !!
    !!
    !! This routine computes \f$ \Delta Q \f$ and updates the solution as \f$ Q^{n+1} = Q^{n} + \Delta Q \f$
    !!
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !------------------------------------------------------------
    type, extends(time_integrator_t), public :: forward_euler_t


    contains
        procedure   :: iterate

        final :: destructor
    end type forward_euler_t
    !-----------------------------------------------------------

contains



    !> Solve for update 'dq'
    !!
    !! \f$ \delta Q = - \delta t R(Q) \f$
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine iterate(self,data,nonlinear_solver,linear_solver,preconditioner)
        class(forward_euler_t),                 intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        class(nonlinear_solver_t),  optional,   intent(inout)   :: nonlinear_solver
        class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner

        character(100)          :: filename
        integer(ik)             :: itime = 1, nsteps, ielem, wcount, iblk, ieqn, idom
        real(rk), allocatable   :: vals(:)


        wcount = 1
        associate ( q => data%sdata%q, dq => data%sdata%dq, rhs => data%sdata%rhs, lhs => data%sdata%lhs, dt => self%dt)

            print*, 'entering time'
            do itime = 1,self%nsteps
                print*, "Step: ", itime


                ! Update Spatial Residual and Linearization (rhs, lin)
                call update_space(data)


                ! Multiply RHS by mass matrix 
                do idom = 1,data%ndomains()
                    do ielem = 1,data%mesh(idom)%nelem

                        do ieqn = 1,data%eqnset(idom)%prop%nprimary_fields()
                            vals = matmul(data%mesh(idom)%elems(ielem)%invmass, rhs%dom(idom)%vecs(ielem)%getvar(ieqn,itime))
                            call rhs%dom(idom)%vecs(ielem)%setvar(ieqn,itime,vals)
                        end do

                    end do !ielem
                end do !idom


                ! Compute update vector
                dq = (-dt) * rhs


                ! Advance solution with update vector
                q  = q + dq



                !
                ! Print diagnostics
                !
                call write_line("   R(Q) - Norm: ", rhs%norm(ChiDG_COMM), delimiter='', io_proc=GLOBAL_MASTER)


                if (wcount == self%nwrite) then
                    write(filename, "(I7,A4)") 1000000+itime, '.plt'
                    call write_tecio_variables_unstructured(data,trim(filename),itime+1)
                    wcount = 0
                end if


                ! Clear residual and linearization storage
                call rhs%clear()
                call lhs%clear()

                wcount = wcount + 1
            end do

        end associate

    end subroutine iterate
    !********************************************************************************************







    
    !>
    !!
    !!
    !!
    !----------------------------------------------------
    subroutine destructor(self)
        type(forward_euler_t),      intent(in) :: self

    end subroutine destructor
    !****************************************************




end module forward_euler










