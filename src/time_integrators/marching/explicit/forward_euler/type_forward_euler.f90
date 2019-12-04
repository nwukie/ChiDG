module type_forward_euler
    use mod_kinds,              only: rk,ik
    use mod_constants,          only: NO_DIFF
    use mod_spatial,            only: update_space
    use type_time_integrator_marching,  only: time_integrator_marching_t
    use type_chidg_data,                only: chidg_data_t
    use type_nonlinear_solver,          only: nonlinear_solver_t
    use type_linear_solver,             only: linear_solver_t
    use type_preconditioner,            only: preconditioner_t
    use type_chidg_vector
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
    !! This routine computes \f$ \Delta Q \f$ and updates the solution 
    !! as \f$ Q^{n+1} = Q^{n} + \Delta Q \f$
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/7/2017
    !!
    !----------------------------------------------------------------------------------------
    type, extends(time_integrator_marching_t), public :: forward_euler_t


    contains

        procedure   :: init
        procedure   :: step

    end type forward_euler_t
    !----------------------------------------------------------------------------------------

contains

    !>  Initialize the forward_euler_t time integrator
    !!
    !!  Create the assembler and atatch it to the time_integrator object so it can 
    !!  be passed to the nonlinear solver
    !!
    !-----------------------------------------------------------------------------------
    subroutine init(self)
        class(forward_euler_t),          intent(inout)   :: self

        call self%set_name('Forward Euler')

    end subroutine init
    !***********************************************************************************






    !> Solve for update 'dq'
    !!
    !! \f$ \delta Q = - \delta t R(Q) \f$
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine step(self,data,nonlinear_solver,linear_solver,preconditioner)
        class(forward_euler_t),                 intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        class(nonlinear_solver_t),  optional,   intent(inout)   :: nonlinear_solver
        class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner

        character(100)          :: filename
        integer(ik)             :: itime, ielem, iblk, ieqn, idom
        real(rk), allocatable   :: vals(:)


        associate ( q   => data%sdata%q,        &
                    dq  => data%sdata%dq,       &
                    rhs => data%sdata%rhs,      &
                    lhs => data%sdata%lhs,      &
                    dt  => data%time_manager%dt)


            !
            ! Update equation Residual (rhs)
            !
            !call update_grid(data)
            call data%update_grid()
            call update_space(data, differentiate=NO_DIFF)
            call self%residual_norm%push_back(rhs%norm(ChiDG_COMM))


            !
            ! Multiply RHS by mass matrix 
            !
            itime = 1
            do idom = 1,data%mesh%ndomains()
                do ielem = 1,data%mesh%domain(idom)%nelem
                    do ieqn = 1,data%eqnset(idom)%prop%nprimary_fields()

                        vals = matmul(data%mesh%domain(idom)%elems(ielem)%invmass, rhs%dom(idom)%vecs(ielem)%getvar(ieqn,itime))
                        call rhs%dom(idom)%vecs(ielem)%setvar(ieqn,itime,vals)

                    end do !ieqn
                end do !ielem
            end do !idom


            !
            ! Compute 'Forward Euler' update vector
            !
            dq = (-dt) * rhs


            !
            ! Advance solution with update vector
            !
            q  = q + dq

!            data%time_manager%t = data%time_manager%t + dt

            !
            ! Clear working vector
            !
            call rhs%clear()


            data%time_manager%t = data%time_manager%t + data%time_manager%dt


        end associate

    end subroutine step
    !*****************************************************************************************




end module type_forward_euler
