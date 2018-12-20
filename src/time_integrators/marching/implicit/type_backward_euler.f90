module type_backward_euler
#include <messenger.h>
    use messenger,                      only: write_line
    use mod_kinds,                      only: rk,ik
    use mod_constants,                  only: ONE, NO_ID
    use mod_spatial,                    only: update_space
!    use mod_update_grid,                only: update_grid

    use type_time_integrator_marching,  only: time_integrator_marching_t
    use type_system_assembler,          only: system_assembler_t

    use type_chidg_data,                only: chidg_data_t
    use type_nonlinear_solver,          only: nonlinear_solver_t
    use type_linear_solver,             only: linear_solver_t
    use type_preconditioner,            only: preconditioner_t
    use type_solver_controller,         only: solver_controller_t
    use type_chidg_vector,              only: chidg_vector_t, sub_chidg_vector_chidg_vector
    use type_chidg_matrix,              only: chidg_matrix_t

    implicit none
    private



    !>  Object implementing the backward Euler time integrator
    !!
    !!  @author Mayank Sharma
    !!  @date   3/31/2017
    !!
    !-------------------------------------------------------------------------------
    type, extends(time_integrator_marching_t),  public  :: backward_euler_t


    contains

        procedure   :: init
        procedure   :: step


    end type backward_euler_t
    !*******************************************************************************



    !>  Object for assembling the implicit system
    !!
    !!  @author Mayank Sharma
    !!  @date   3/31/2017
    !!
    !-------------------------------------------------------------------------------
    type, extends(system_assembler_t),  public  :: assemble_backward_euler_t

        type(chidg_vector_t)    :: q_n

    contains

        procedure   :: assemble


    end type assemble_backward_euler_t
    !*******************************************************************************




    !>  Control the lhs update inside the nonlinear solver.
    !!
    !!  Reference:
    !!  Persson, P.-O., "High-Order Navier-Stokes Simulations using a Sparse 
    !!  Line-Based Discontinuous Galerkin Method", AIAA-2012-0456
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/22/2017
    !!
    !--------------------------------------------------------------------------------
    type, extends(solver_controller_t), public :: backwardeuler_solver_controller_t

    contains

        procedure   :: update_lhs

    end type backwardeuler_solver_controller_t
    !********************************************************************************
















contains



    !>  Initialize the backward_euler_t time integrator
    !!
    !!  Create the assembler and attach it to the time_integrator object so it can
    !!  be passed to the nonlinear solver
    !!
    !!  @author Mayank Sharma
    !!  @date   3/31/2017
    !!
    !-------------------------------------------------------------------------------
    subroutine init(self)
        class(backward_euler_t),    intent(inout)   :: self

        integer(ik)                     :: ierr
        type(assemble_backward_euler_t) :: assemble_backward_euler

        call self%set_name('Backward Euler')

        if (allocated(self%system)) deallocate(self%system)
        allocate(self%system, source=assemble_backward_euler, stat=ierr)
        if (ierr /= 0) call AllocationError


    end subroutine init
    !*******************************************************************************



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
    !--------------------------------------------------------------------------------------
    subroutine step(self,data,nonlinear_solver,linear_solver,preconditioner)
        class(backward_euler_t),                intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        class(nonlinear_solver_t),  optional,   intent(inout)   :: nonlinear_solver
        class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner

        class(*),   allocatable     :: assemble_type

        type(backwardeuler_solver_controller_t),    save    :: solver_controller


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
        !call update_grid(data)
        call data%update_grid()
        call nonlinear_solver%solve(data,self%system,linear_solver,preconditioner)

        !
        ! Store end residual from nonlinear solver
        !
        call self%residual_norm%push_back(nonlinear_solver%residual_norm%at(nonlinear_solver%residual_norm%size()))


    end subroutine step
    !**************************************************************************************



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
    !--------------------------------------------------------------------------------------
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


        associate ( q   => data%sdata%q,   &
                    dq  => data%sdata%dq,  &
                    lhs => data%sdata%lhs, &
                    rhs => data%sdata%rhs)

        !
        ! Clear data containers
        !
        call rhs%clear()
        if (differentiate) call lhs%clear()


        !
        ! Get spatial update
        !
        call update_space(data,differentiate,timing)


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
    !**************************************************************************************







    !>  Control algorithm for selectively updating the lhs matrix in the
    !!  nonlinear solver.
    !!
    !!  Reference:
    !!  Persson, P.-O., "High-Order Navier-Stokes Simulations using a Sparse 
    !!  Line-Based Discontinuous Galerkin Method", AIAA-2012-0456
    !!
    !!  @author Nathan A. Wukie
    !!  @date   6/22/2017
    !!
    !!  @param[in]  niter               Number of newton iterations
    !!  @param[in]  residual_ratio      R_{i}/R_{i-1}
    !!
    !----------------------------------------------------------------------------
    function update_lhs(self,A,niter,residual_ratio) result(update)
        class(backwardeuler_solver_controller_t),   intent(inout)   :: self
        type(chidg_matrix_t),                       intent(in)      :: A
        integer(ik),                                intent(in)      :: niter
        real(rk),                                   intent(in)      :: residual_ratio

        logical :: update

        ! Update lhs if:
        !   1: If matrix(lhs/A) hasn't been updated before
        !   2: number of newton iterations > 10
        !   3: residual norm increases by factor of 10 (divergence)
        !   4: being forced
        if ( all(A%stamp == NO_ID)      .or. &
            (niter > 10)                .or. &
            (residual_ratio > 10._rk)   .or. &
            (self%force_update_lhs) ) then
            update = .true.
        else
            update = .false.
        end if

        ! Store action
        self%lhs_updated = update

        ! Turn off forced update
        self%force_update_lhs = .false.

    end function update_lhs
    !****************************************************************************


















end module type_backward_euler
