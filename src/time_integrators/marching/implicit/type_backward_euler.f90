module type_backward_euler
#include <messenger.h>
    use messenger,                      only: write_line
    use mod_kinds,                      only: rk,ik
    use mod_constants,                  only: ONE
    use mod_spatial,                    only: update_space
    use mod_update_grid,                only: update_grid

    use type_time_integrator_marching,  only: time_integrator_marching_t
    use type_system_assembler,          only: system_assembler_t

    use type_chidg_data,                only: chidg_data_t
    use type_nonlinear_solver,          only: nonlinear_solver_t
    use type_linear_solver,             only: linear_solver_t
    use type_preconditioner,            only: preconditioner_t
    use type_chidg_vector,              only: chidg_vector_t, sub_chidg_vector_chidg_vector

    implicit none
    private



    !>  Object implementing the backward Euler time integrator
    !!
    !!  @author Mayank Sharma
    !!  @date   3/31/2017
    !!
    !------------------------------------------------------------------------------------------------
    type, extends(time_integrator_marching_t),  public  :: backward_euler_t


    contains

        procedure   :: init
        procedure   :: step


    end type backward_euler_t
    !************************************************************************************************



    !>  Object for assembling the implicit system
    !!
    !!  @author Mayank Sharma
    !!  @date   3/31/2017
    !!
    !------------------------------------------------------------------------------------------------
    type, extends(system_assembler_t),  public  :: assemble_backward_euler_t

        type(chidg_vector_t)    :: q_n

    contains

        procedure   :: assemble


    end type assemble_backward_euler_t
    !************************************************************************************************



contains



    !>  Initialize the backward_euler_t time integrator
    !!
    !!  Create the assembler and attach it to the time_integrator object so it can
    !!  be passed to the nonlinear solver
    !!
    !!  @author Mayank Sharma
    !!  @date   3/31/2017
    !!
    !------------------------------------------------------------------------------------------------
    subroutine init(self,data)
        class(backward_euler_t),    intent(inout)   :: self
        type(chidg_data_t),         intent(in)      :: data

        integer(ik)                     :: ierr
        type(assemble_backward_euler_t) :: assemble_backward_euler


        if (allocated(self%system)) deallocate(self%system)
        allocate(self%system, source=assemble_backward_euler, stat=ierr)
        if (ierr /= 0) call AllocationError


    end subroutine init
    !************************************************************************************************



    !>  Solution advancement via the backward Euler method
    !!
    !!  Given the system of partial differential equations consisting of the time-derivative of the
    !!  the solution vector and a spatial residual as
    !!
    !!  \f$ M \frac{\partial Q}{\partial t} + R(Q) = 0 \f$
    !!
    !!  the time derivative is discretized with a bacward difference 
    !!
    !!  \f$ M \frac{Q^{n + 1} - Q^{n}}{\Delta t} + R(Q^{n + 1}) = 0 \f$
    !!
    !!  which yields a nonlinear equation
    !!
    !!  \f$  \frac{\Delta Q}{\Delta t}M + R(Q^{n + 1}) = 0 \f$
    !!
    !!  which results in the Newton iteration
    !!
    !!  \f$ \left(\frac{M}{\Delta t} + \frac{\partial R(Q^{m})}{\partial Q}\right) \delta Q^{m} = 
    !!      -\frac{M}{\Delta t}\Delta Q^{m} - R(Q^{m}) \f$
    !!  \f$ Q^{m} = Q^{n} + \Delta Q^{m} \f$
    !!
    !!  where \f$ \delta Q^{m} = \Delta Q^{m + 1} -  \Delta Q^{m} \f$ for the mth Newton iteration
    !!
    !!  @author Mayank Sharma
    !!  @date   3/31/2017
    !!
    !------------------------------------------------------------------------------------------------
    subroutine step(self,data,nonlinear_solver,linear_solver,preconditioner)
        class(backward_euler_t),                intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        class(nonlinear_solver_t),  optional,   intent(inout)   :: nonlinear_solver
        class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner

        class(*),   allocatable     :: assemble_type


        select type(associate_name => self%system)
            type is (assemble_backward_euler_t)

                !
                ! Store solution at nth time step to a separate vector
                !
                associate_name%q_n = data%sdata%q
                

        end select


        !
        ! Solve assembled nonlinear system, the nonlinear update is the step in time
        ! System assembled in subroutine assemble
        !
        data%time_manager%t = data%time_manager%t + data%time_manager%dt
        call update_grid(data)
        call nonlinear_solver%solve(data,self%system,linear_solver,preconditioner)

        !
        ! Store end residual from nonlinear solver
        !
        call self%residual_norm%push_back(nonlinear_solver%residual_norm%at(nonlinear_solver%residual_norm%size()))


    end subroutine step
    !************************************************************************************************



    !>  Assemble the system for the backward Euler equations with temporal contributions
    !!
    !!  @author Mayank Sharma
    !!  @date   3/31/2017
    !!
    !!  \f$ lhs = \frac{\partial R(Q)}{\partial Q} \f$
    !!  \f$ rhs = R(Q) \f$
    !!  \f$ M   = element mass matrix \f$
    !!
    !!  For system assembly with temporal contributions
    !!
    !!  \f$ lhs = \frac{M}{dt} + lhs \f$
    !!  \f$ rhs = \frac{M}{dt}*dq + rhs \f$
    !!
    !------------------------------------------------------------------------------------------------
    subroutine assemble(self,data,differentiate,timing)
        class(assemble_backward_euler_t),   intent(inout)               :: self
        type(chidg_data_t),                 intent(inout)               :: data
        logical,                            intent(in)                  :: differentiate
        real(rk),                           intent(inout),  optional    :: timing

        integer(ik)                 :: ntime, itime, idom, ielem, ivar, imat, ierr, &
                                       nterms, rstart, rend, cstart, cend
        real(rk)                    :: dt
        type(chidg_vector_t)        :: delta_q
        real(rk),   allocatable     :: temp_1(:), temp_2(:)


        !
        ! Get spatial update
        !
        call update_space(data,differentiate,timing)


        associate ( q   => data%sdata%q,   &
                    dq  => data%sdata%dq,  &
                    lhs => data%sdata%lhs, &
                    rhs => data%sdata%rhs)

            !
            ! Get no. of time levels ( = 1 for time marching) and time step
            !
            ntime = data%time_manager%ntime
            dt    = data%time_manager%dt

            
            !
            ! Compute \f$ q^{m} - q^{n} \f$
            ! Used to assemble rhs
            !
            delta_q = sub_chidg_vector_chidg_vector(q,self%q_n)


            do itime = 1,ntime
                do idom = 1,data%mesh%ndomains()

                    !
                    ! Allocate temporary arrays 
                    !
                    if (allocated(temp_1) .and. allocated(temp_2)) deallocate(temp_1,temp_2)
                    allocate(temp_1(data%mesh%domain(idom)%nterms_s), temp_2(data%mesh%domain(idom)%nterms_s), stat=ierr)
                    if (ierr /= 0) call AllocationError

                    do ielem = 1,data%mesh%domain(idom)%nelem
                        do ivar = 1,data%eqnset(idom)%prop%nprimary_fields()

                            !
                            ! Assemble lhs 
                            !
                            nterms = data%mesh%domain(idom)%elems(ielem)%nterms_s
                            rstart = 1 + (ivar - 1)*nterms
                            rend   = (rstart - 1) + nterms
                            cstart = rstart
                            cend   = rend

                            ! Add mass matrix divided by dt to the block diagonal
                            imat   = lhs%dom(idom)%lblks(ielem,itime)%get_diagonal()
                            lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(rstart:rend,cstart:cend) = lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(rstart:rend,cstart:cend) + &
                            data%mesh%domain(idom)%elems(ielem)%mass*(ONE/dt)

    
                            !
                            ! Assemble rhs
                            !
                            temp_1 = (ONE/dt)*matmul(data%mesh%domain(idom)%elems(ielem)%mass, delta_q%dom(idom)%vecs(ielem)%getvar(ivar,itime))
                            temp_2 = rhs%dom(idom)%vecs(ielem)%getvar(ivar,itime) + temp_1
                            call rhs%dom(idom)%vecs(ielem)%setvar(ivar,itime,temp_2)
                                   
                        end do  ! ivar
                    end do  ! ielem

                end do  ! idom
            end do  ! itime

        end associate


    end subroutine assemble
    !************************************************************************************************




















end module type_backward_euler
