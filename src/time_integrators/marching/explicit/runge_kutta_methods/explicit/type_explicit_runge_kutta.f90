module type_explicit_runge_kutta
    use mod_kinds,                      only: rk,ik
    use mod_constants,                  only: ZERO,HALF,ONE,TWO,SIX
    use mod_spatial,                    only: update_space

    use type_time_integrator_marching,  only: time_integrator_marching_t
    use type_chidg_data,                only: chidg_data_t
    use type_nonlinear_solver,          only: nonlinear_solver_t
    use type_linear_solver,             only: linear_solver_t
    use type_preconditioner,            only: preconditioner_t
    use type_chidg_vector,           
    
    use mod_define_explicit_RK_methods, only: method_selector
    
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
    type, extends(time_integrator_marching_t), public :: explicit_runge_kutta_t



    contains

        procedure   :: init
        procedure   :: step


    end type explicit_runge_kutta_t
    !*********************************************************************************************************************



contains



    !>  Initialize the forward_euler_t time integrator
    !!
    !!  Create the assembler and atatch it to the time_integrator object so it can 
    !!  be passed to the nonlinear solver
    !!
    !-----------------------------------------------------------------------------------
    subroutine init(self)
        class(explicit_runge_kutta_t),          intent(inout)   :: self

        call self%set_name('ERK')

    end subroutine init
    !***********************************************************************************



    !>  Solve for update dq
    !!
    !!  \f$ \Delta Q_{i} = -\Delta t m^{-1} R(Q^{n} + \sum_{j = 1}^{i-1} a_{ij} \Delta Q_{j}) \f$
    !!
    !!  @author Mayank Sharma
    !!  @date   18/1/2017
    !!
    !!
    !---------------------------------------------------------------------------------------------------------------------
    subroutine step(self,data,nonlinear_solver,linear_solver,preconditioner)
        class(explicit_runge_kutta_t),              intent(inout)   :: self
        type(chidg_data_t),                         intent(inout)   :: data
        class(nonlinear_solver_t),      optional,   intent(inout)   :: nonlinear_solver
        class(linear_solver_t),         optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),        optional,   intent(inout)   :: preconditioner

        !
        ! Variables read from subroutines defining the explicit RK methods
        !
        integer(ik)                          :: nstage   
        real(rk),            allocatable     :: a(:,:)   
        real(rk),            allocatable     :: b(:)
        !
        ! Variables defined for use in this module
        !
        type(chidg_vector_t), allocatable    :: delq(:)
        type(chidg_vector_t)                 :: q_n, q_adv
        character(100)                       :: filename
        character(len = :),  allocatable     :: time_scheme
        integer(ik)                          :: itime = 1, nsteps, ielem, wcount, iblk, ieqn, idom, istage, j
        real(rk),            allocatable     :: vals(:)


        associate( q   => data%sdata%q,      &
                   dq  => data%sdata%dq,     &
                   rhs => data%sdata%rhs,    &
                   lhs => data%sdata%lhs,    & 
                   dt  => data%time_manager%dt )
        
            !
            ! Get name of the RK method being used
            !
            time_scheme = data%time_manager%get_name()

            !
            ! Get nstage, a and b for a particular scheme
            ! All methods defined in mod_explicit_RK_methods.f90 
            !
            call method_selector(time_scheme,nstage,a,b)


            !
            ! Copy the solution vector at the present time level to a temporary vector
            !
            q_n = q

            !
            ! q_adv is a temporary vector used to advance the solution after each stage
            ! Computed using an accumulated sum
            ! Initialize q_adv
            !
            q_adv = q_n

            if (allocated(delq)) deallocate(delq)
            allocate(delq(nstage))

            !
            ! Loop through the stages of the present RK method
            !
            do istage = 1,nstage
                
                !
                ! Update Spatial Residual and Linearization (rhs, lin)
                !
                call update_space(data,differentiate=.false.)
                call self%residual_norm%push_back(rhs%norm(ChiDG_COMM))

                do idom = 1,data%mesh%ndomains()
                    do ielem = 1,data%mesh%domain(idom)%nelem
                        do ieqn = 1,data%eqnset(idom)%prop%nprimary_fields()
                            
                            vals = matmul(data%mesh%domain(idom)%elems(ielem)%invmass, &
                                          rhs%dom(idom)%vecs(ielem)%getvar(ieqn,itime))
                            call rhs%dom(idom)%vecs(ielem)%setvar(ieqn,itime,vals)

                        end do  ! ieqn
                    end do  ! ielem
                end do  ! idom

                !
                ! Calculate stagewise update and its contribution to the final update
                !
                dq = (-dt)*rhs
                delq(istage) = b(istage)*dq

                !
                ! Compute the solution vector to pass to update_space for the next stage
                ! Refer to scheme description above
                !
                if (istage < nstage) then
                    do j = 1,istage
                        q = q + (a(istage + 1,j)*delq(j)/b(j))
                    end do
                end if

                !
                ! Clear working vector
                !
                call rhs%clear()

                !
                ! Update q_adv after each stage
                !
                q_adv = add_chidg_vector_chidg_vector(q_adv, delq(istage))

            end do  ! istage


            !
            ! Advance solution
            ! 
            q = q_adv


        end associate


    end subroutine step
    !*********************************************************************************************************************




















end module type_explicit_runge_kutta
