module type_steady
#include <messenger.h>
    use mod_kinds,                      only: rk,ik
    use mod_constants,                  only: ZERO, dQ_DIFF, NO_DIFF
    use mod_spatial,                    only: update_space
    use mod_update_functionals,         only: update_functionals
    use mod_chidg_mpi,                  only: GLOBAL_MASTER

    use type_time_integrator_steady,    only: time_integrator_steady_t
    use type_system_assembler,          only: system_assembler_t

    use type_chidg_data,                only: chidg_data_t
    use type_nonlinear_solver,          only: nonlinear_solver_t
    use type_linear_solver,             only: linear_solver_t
    use type_preconditioner,            only: preconditioner_t
    use type_chidg_vector
    implicit none
    private



    !>  Object implementing solving a system of equations to steady-state.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/29/2016
    !!
    !----------------------------------------------------------------------------------------
    type, extends(time_integrator_steady_t), public :: steady_t

    contains
        procedure   :: init
        procedure   :: step
        procedure   :: compute_adjoint
        procedure   :: compute_functionals
    end type steady_t
    !****************************************************************************************




    !>  Object for assembling the implicit system.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2017
    !!
    !----------------------------------------------------------------------------------------
    type, extends(system_assembler_t), public :: assemble_steady_t

    contains
        procedure   :: assemble
    end type assemble_steady_t
    !****************************************************************************************



contains





    !>  Initialize the steady_t time integrator.
    !!
    !!  Main activity here, create the assembler and attach it to the time_integrator 
    !!  object so it can be passed to the nonlinear solver.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2017
    !!
    !----------------------------------------------------------------------------------------
    subroutine init(self)
        class(steady_t),    intent(inout)   :: self

        integer(ik)             :: ierr
        type(assemble_steady_t) :: assemble_steady

        ! Set name
        call self%set_name('Steady')

        ! Allocate assembler
        if (allocated(self%system)) deallocate(self%system)
        allocate(self%system, source=assemble_steady, stat=ierr)
        if (ierr /= 0) call AllocationError

    end subroutine init
    !*****************************************************************************************




    !>  Steady equation set. No time.
    !!
    !!  Solving R(Q) = 0
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/29/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine step(self,data,nonlinear_solver,linear_solver,preconditioner)
        class(steady_t),                        intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        class(nonlinear_solver_t),  optional,   intent(inout)   :: nonlinear_solver
        class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner


        ! Simply solve the nonlinear system. No iteration in time.
        call nonlinear_solver%solve(data,self%system,linear_solver,preconditioner)


        ! Store end residual from nonlinear solver.
        call self%residual_norm%push_back(nonlinear_solver%residual_norm%at(nonlinear_solver%residual_norm%size()))


    end subroutine step
    !******************************************************************************************





    !>  Assemble the system for the 'Steady' equations with no time derivative contributions.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine assemble(self,data,differentiate,timing)
        class(assemble_steady_t),   intent(inout)               :: self
        type(chidg_data_t),         intent(inout)               :: data
        integer(ik),                intent(in)                  :: differentiate
        real(rk),                   intent(inout),  optional    :: timing

        integer(ik) :: eqn_ID, idom, ielem, itime, imat, ifield, irow_start, icol_start, nterms, nfields, ifield_row, ifield_col
        real(rk),   allocatable     :: elem_field(:), elem_res(:)

        call data%sdata%rhs%clear()
        if (differentiate == dQ_DIFF) call data%sdata%lhs%clear()

        ! Steady equation, so we only need the spatial operators computed.
        data%time_manager%itime = 1
        call data%update_grid()
        !call update_grid(data)

        call update_space(data,differentiate,timing)

    end subroutine assemble
    !******************************************************************************************





    !>  Computation of adjoint variables for steady case (1 chidg_vector of adjoint variables).
    !!
    !!  @author Matteo Ugolotti
    !!  @date   17/6/2017
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute_adjoint(self,data,linear_solver,preconditioner)
        class(steady_t),                        intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner
        
        integer(ik)                 :: istep, ifunc
        character(:),   allocatable :: func_name

        ! istep = 1 for steady
        istep = data%time_manager%istep

        associate( lhs  => data%sdata%lhs,               &
                   adj  => data%sdata%adjoint%v,         &
                   Jq   => data%sdata%adjoint%Jq,        &
                   adjoint => data%sdata%adjoint )

            ! Copy the unique solution vector q_time in q to carry out the update_spatial
            call data%sdata%q%assemble()
            associate(q => data%sdata%q, q_time => data%sdata%adjoint%q_time(istep))
                q = q_time
                !data%sdata%q = data%sdata%adjoint%q_time(istep)
            end associate

            ! Assemble the system updating spatial residuals and linearization (rhs,lhs)
            call self%system%assemble(data,differentiate=dQ_DIFF)

            ! Update functional/s
            call update_functionals(data,differentiate=dQ_DIFF)

            ! Loop through functionals and compute adjoint for each of them
            do ifunc = 1,size(Jq)

                ! Starting adjoint computation message for ifunc
                func_name = data%functional_group%fcl_entities(ifunc)%func%get_name()
                call write_line('-  Computing primary adjoint variables for',func_name, 'functional...',io_proc=GLOBAL_MASTER)

                ! Compute adjoint variables for ifunc and store linear solver info
                call linear_solver%solve(lhs,adj(ifunc,istep),Jq(ifunc),preconditioner)

                call adjoint%store_solver_info(ifunc,istep,linear_solver%niter,             &
                                                           linear_solver%timer%elapsed(),   &
                                                           linear_solver%residual_norm)

                ! Done adjoint computation message for ifunc
                call write_line('-  Primary adjoint variables for',func_name, 'functional computed.',io_proc=GLOBAL_MASTER)

            end do

        end associate

    end subroutine compute_adjoint
    !******************************************************************************************





    !>  Computation of functionals
    !!
    !!  @author Matteo Ugolotti
    !!  @date   6/17/2017
    !!
    !!  TODO: maybe move this to time_integrator_t, it is the same for each time integrator
    !!        steady or unsteady since it is based on the a single solution in time. 
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute_functionals(self,data)
        class(steady_t),                        intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        
        ! Assemble the system updating spatial residuals and linearization (rhs,lhs)
        !
        ! We might not need this again
        ! call self%system%assemble(data,differentiate=dQ_DIFF)

        ! Update functional
        call update_functionals(data,differentiate=NO_DIFF)

    end subroutine compute_functionals
    !******************************************************************************************



    
end module type_steady
