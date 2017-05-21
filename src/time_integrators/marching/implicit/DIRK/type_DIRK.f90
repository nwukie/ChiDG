module type_DIRK
#include<messenger.h>
    use messenger,                      only: write_line
    use mod_kinds,                      only: rk, ik
    use mod_constants,                  only: ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX, SEVEN
    use mod_spatial,                    only: update_space

    use type_time_integrator_marching,  only: time_integrator_marching_t
    use type_system_assembler,          only: system_assembler_t

    use type_chidg_data,                only: chidg_data_t
    use type_nonlinear_solver,          only: nonlinear_solver_t           
    use type_linear_solver,             only: linear_solver_t
    use type_preconditioner,            only: preconditioner_t
    use type_chidg_vector

    implicit none
    private


    !>  Object implementing the diagonally implicit RK time integrator
    !!
    !!  @author Mayank Sharma
    !!  @date   5/20/2017
    !!
    !-----------------------------------------------------------------------------------------------------------------
    type, extends(time_integrator_marching_t),  public      :: DIRK_t


    contains

        procedure   :: init
        procedure   :: step


    end type DIRK_t
    !*****************************************************************************************************************




    !>  Object for assembling the implicit system
    !!
    !!  @author Mayank Sharma
    !!  @date   5/20/2017
    !!
    !-----------------------------------------------------------------------------------------------------------------
    type, extends(system_assembler_t),  public      :: assemble_DIRK_t

        type(chidg_vector_t)    :: q_n
        type(chidg_vector_t)    :: q_n_stage

    contains

        procedure   :: assemble


    end type assemble_DIRK_t
    !*****************************************************************************************************************



contains




    !>  Initialize the DIRK_t time integrator
    !!
    !!  Create the assembler and atatch it to the time_integrator object so it can 
    !!  be passed to the nonlinear solver
    !!
    !!  @author Mayank Sharma
    !!  @date   5/20/2017
    !!
    !-----------------------------------------------------------------------------------------------------------------
    subroutine init(self,data)
        class(DIRK_t),          intent(inout)   :: self
        type(chidg_data_t),     intent(in)      :: data

        integer(ik)             :: ierr
        type(assemble_DIRK_t)   :: assemble_DIRK


        if (allocated(self%system)) deallocate(self%system)
        allocate(self%system, source=assemble_DIRK, stat=ierr)
        if (ierr /= 0) call AllocationError


    end subroutine init
    !*****************************************************************************************************************




    !>
    !!
    !!  @author Mayank Sharma
    !!  @date   5/20/2017
    !!
    !-----------------------------------------------------------------------------------------------------------------
    subroutine step(self,data,nonlinear_solver,linear_solver,preconditioner)
        class(DIRK_t),                          intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        class(nonlinear_solver_t),  optional,   intent(inout)   :: nonlinear_solver
        class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner

        integer(ik),    parameter   :: nstage = 3
        type(chidg_vector_t)        :: dq(nstage), q_temp, q_n
        real(rk)                    :: A(nstage,nstage), b(nstage)
        real(rk)                    :: alpha, tau, b1, b2
        integer(ik)                 :: istage, jstage


        !
        ! Constants used to define DIRK coefficients
        !
        alpha = 0.435866521508459
        tau   = (ONE + alpha)/TWO
        b1    = -(SIX*(alpha*alpha) - (16._rk*alpha) + ONE)/FOUR
        b2    = (SIX*(alpha*alpha) - (20._rk*alpha) + FIVE)/FOUR


        !
        ! DIRK coefficient arrays
        !
        A     = transpose(reshape([alpha, ZERO,  ZERO, &
                                   tau,   alpha, ZERO, &
                                   b1,    b2,    alpha], shape(A))) 
        b     = [b1, b2, alpha]

        !
        ! Store solution at nth time step to a separate vector for use in this subroutine
        !
        q_n = data%sdata%q

        select type(associate_name => self%system)
            type is (assemble_DIRK_t)
                associate_name%q_n = data%sdata%q

        end select


        !
        ! Initialize update vector array
        !
        do istage = 1, nstage

            call dq(istage)%init(data%mesh, data%time_manager%ntime)
            call dq(istage)%set_ntime(data%time_manager%ntime)
            call dq(istage)%clear()

        end do


        associate ( q => data%sdata%q )
            
            do istage = 1, nstage

                select case(istage)
                    case(1)
                        select type(an => self%system)
                            type is (assemble_DIRK_t)
                                call an%q_n_stage%init(data%mesh,data%time_manager%ntime)
                                call an%q_n_stage%set_ntime(data%time_manager%ntime)
                                call an%q_n_stage%clear()

                        end select
                        q_temp = q_n
                    case(2)
                        select type(an => self%system)
                            type is (assemble_DIRK_t)
                                an%q_n_stage = (tau - alpha)*dq(1)

                        end select
                        q_temp = q_n + (tau - alpha)*dq(1)
                    case(3)
                        select type(an => self%system)
                            type is (assemble_DIRK_t)
                                an%q_n_stage = b1*dq(1) + b2*dq(2)

                        end select
                        q_temp = q_n + b1*dq(1) + b2*dq(2)

                end select

                !
                ! Solve assembled nonlinear system, the nonlinear update is the stagewise update
                ! System assembled in subroutine assemble
                !
                call nonlinear_solver%solve(data,self%system,linear_solver,preconditioner)

                !
                ! Store end residual from nonlinear solver
                !
                call self%residual_norm%push_back(nonlinear_solver%residual_norm%at(nonlinear_solver%residual_norm%size()))

                !
                ! Store stagewise update
                !
                dq(istage) = (q - q_temp)/alpha

            end do

        end associate


    end subroutine step
    !*****************************************************************************************************************




    !>  Assemble the system for the DIRK equations with temporal contributons
    !!
    !!  @author Mayank Sharma
    !!  @date   5/20/2017
    !!
    !-----------------------------------------------------------------------------------------------------------------
    subroutine assemble(self,data,timing,differentiate)
        class(assemble_DIRK_t), intent(inout)               :: self
        type(chidg_data_t),     intent(inout)               :: data
        real(rk),               intent(inout),  optional    :: timing
        logical,                intent(in),     optional    :: differentiate

        real(rk),   parameter       :: alpha = 0.435866521508459
        type(chidg_vector_t)        :: delta_q 
        real(rk)                    :: dt
        integer(ik)                 :: ntime, itime, idom, ielem, ivar, imat, ierr, &       
                                       nterms, rstart, rend, cstart, cend
        real(rk),   allocatable     :: temp_1(:), temp_2(:)


        associate( q   => data%sdata%q,   &
                   dq  => data%sdata%dq,  &
                   lhs => data%sdata%lhs, & 
                   rhs => data%sdata%rhs)

            !
            ! Get spatial update
            !
            call update_space(data,timing,differentiate)

            !
            ! Get no. of time levels (=1 for time marching) and time step
            !
            ntime = data%time_manager%ntime
            dt    = data%time_manager%dt

            !
            ! Compute \f$ \Delta Q^{m}_{i}\f$
            ! Used to assemble rhs
            !
            delta_q = (q - self%q_n - self%q_n_stage)/alpha

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
                            rend   = (rstart - 1) +  nterms
                            cstart = rstart
                            cend   = rend

                            ! Add mass matrix divided by (alpha*dt) to the block diagonal
                            imat = lhs%dom(idom)%lblks(ielem,itime)%get_diagonal()
                            lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(rstart:rend,cstart:cend) = lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(rstart:rend,cstart:cend) + &
                            data%mesh%domain(idom)%elems(ielem)%mass*(ONE/(alpha*dt))


                            !
                            ! Assemble rhs
                            !
                            temp_1 = (ONE/(alpha*dt))*matmul(data%mesh%domain(idom)%elems(ielem)%mass, delta_q%dom(idom)%vecs(ielem)%getvar(ivar,itime))
                            temp_2 = rhs%dom(idom)%vecs(ielem)%getvar(ivar,itime) + temp_1
                            call rhs%dom(idom)%vecs(ielem)%setvar(ivar,itime,temp_2)

                        end do  ! ivar
                    end do  ! ielem

                end do  ! idom
            end do  ! itime

        end associate


    end subroutine assemble
    !*****************************************************************************************************************




















end module type_DIRK
