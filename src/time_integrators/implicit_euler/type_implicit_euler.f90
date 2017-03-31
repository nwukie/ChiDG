module type_implicit_euler
#include <messenger.h>
    use messenger,                      only: write_line
    use mod_kinds,                      only: rk,ik
    use mod_constants,                  only: ONE
    use mod_spatial,                    only: update_space

    use type_time_integrator_marching,  only: time_integrator_marching_t
    use type_system_assembler,          only: system_assembler_t

    use type_chidg_data,                only: chidg_data_t
    use type_nonlinear_solver,          only: nonlinear_solver_t
    use type_linear_solver,             only: linear_solver_t
    use type_preconditioner,            only: preconditioner_t
    use type_chidg_vector,              only: chidg_vector_t

    implicit none
    private



    !>  Object implementing the backward Euler time integrator
    !!
    !!  @author Mayank Sharma
    !!  @date   3/31/2017
    !!
    !------------------------------------------------------------------------------------------------
    type, extends(time_integrator_marching_t),  public  :: implicit_euler_t


    contains

        procedure   :: init
        procedure   :: step


    end type implciit_euler_t
    !************************************************************************************************



    !>  Object for assembling the implicit system
    !!
    !!  @author Mayank Sharma
    !!  @date   3/31/2017
    !!
    !------------------------------------------------------------------------------------------------
    type, extends(system_assembler_t),  public  :: assemble_implicit_euler_t
    

    contains

        procedure   :: assemble


    end type assemble_implicit_euler_t
    !************************************************************************************************



contains



    !>  Initialize the implicit_euler_t time integrator
    !!
    !!  Create the assembler and attach it to the time_integrator object so it can
    !!  be passed to the nonlinear solver
    !!
    !!  @author Mayank Sharma
    !!  @date   3/31/2017
    !!
    !------------------------------------------------------------------------------------------------
    subroutine init(self,data)
        class(implicit_euler_t),    intent(inout)   :: self
        type(chidg_data_t),         intent(in)      :: data

        integer(ik)                     :: ierr
        type(assemble_implicit_euler_t) :: assemble_implicit_euler


        if (allocated(self%system)) deallocate(self%system)
        allocate(self%system, source=assemble_implicit_euler, stat=ierr)
        if (ierr /= 0) call AllocationError


    end subroutine init
    !************************************************************************************************



    !>
    !!
    !!  @author Mayank Sharma
    !!  @date   3/31/2017
    !!
    !------------------------------------------------------------------------------------------------
    subroutine step(self,data,nonlinear_solver,linear_solver,preconditioner)
        class(implicit_euler_t),                intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        class(nonlinear_solver_t),  optional,   intent(inout)   :: nonlinear_solver
        class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner


    end subroutine step
    !************************************************************************************************



    !>  Assemble the system for the implicit Euler equations with temporal contributions
    !!
    !!  @author Mayank Sharma
    !!  @date   3/31/2017
    !!
    !------------------------------------------------------------------------------------------------
    subroutine assemble(self,data,timing,differentiate)
        class(assemble_implicit_euler_t),   intent(inout)               :: self
        type(chidg_data_t),                 intent(inout)               :: data
        real(rk),                           intent(inout),  optional    :: timing
        logical,                            intent(in),     optional    :: differentiate

        integer(ik)     :: ntime, itime, idom, ielem, ivar, imat, ierr, &
                           nterms, rstart, rend, cstart, cend
        real(rk)        :: dt


        !
        ! Get spatial update
        !
        call update_space(data,timing,differentiate)

        associate ( q   => data%sdata%q,   &
                    dq  => data%sdata%dq,  &
                    lhs => data%sdata%lhs, &
                    rhs => data%sdata%rhs)

            ntime = data%time_manager%ntime
            dt    = data%time_manager%dt

            do itime = 1,ntime
                do idom = 1,data%ndomains()

                    do ielem = 1,data%mesh(idom)%nelem
                        do ivar = 1,data%eqnset(idom)%prop%nprimary_fields()

                            nterms = data%mesh(idom)%elems(ielem)%nterms_s
                            rstart = 1 + (ieqn - 1)*nterms
                            rend   = (rstart - 1) + nterms
                            cstart = rstart
                            cend   = rend

                            ! Add mass matrix divided by dt to the block diagonal
                            imat   = lhs%dom(idom)%lblks(ielem,itime)%get_diagonal()
                            lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(rstart:rend,cstart:cend) = &
                            lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(rstart:rend,cstart:cend) + &
                            data%mesh(idom)%elems(ielem)%mass*(ONE/dt)
                                   
                        end do  ! ivar
                    end do  ! ielem

                end do  ! idom
            end do  ! itime

        end associate


    end subroutine assemble
    !************************************************************************************************




















end module type_implicit_euler
