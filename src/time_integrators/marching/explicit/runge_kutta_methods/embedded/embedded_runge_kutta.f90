module embedded_runge_kutta
    use mod_kinds,                          only: rk,ik
    use mod_constants,                      only: ZERO,TENTH,HALF,ONE,FOUR
    use type_time_integrator,               only: time_integrator_t
    use type_chidg_data,                    only: chidg_data_t
    use type_nonlinear_solver,              only: nonlinear_solver_t 
    use type_linear_solver,                 only: linear_solver_t
    use type_preconditioner,                only: preconditioner_t
    use type_chidgVector,

    use mod_spatial,                        only: update_space
 
    use mod_tecio,                          only: write_tecio_variables_unstructured

    use mod_define_embedded_RK_methods,     only: method_selector

    implicit none
    private


    !>  Solution advancement via a general embedded Runge-Kutta method
    !!
    !!  The system of partial differential equations consisting of the time-derivative of the
    !!  solution vector and a spatial residual is given as
    !!
    !!  \f$ \frac{\partial Q}{\partial t} + R(Q) = 0 \f$
    !!
    !!  The solution of the next time level is then computed as
    !!
    !!  \f$ Q^{n+1} = Q^{n} + \Delta Q \f$
    !!
    !!  where
    !!
    !!  \f$ \Delta Q = \sum_{i = 1}^{s} b_{i} \Delta Q_{i} \f$
    !!
    !!  with 
    !!    
    !!  \f$ \Delta Q_{i} = -\Delta t m^{-1} R(Q^{n} + \sum_{j = 1}^{i - 1} a_{ij} \Delta Q_{j}) \f$
    !!
    !!  where \f$ M \f$ refers to the mass matrix, \f$ s \f$ refers to the number of stages in the
    !!  method and \f$ a_{ij}'s & b_{j}'s \f$ refer to the coefficients defined by the Butcher tableau
    !!
    !!  In an embedded RK method, the solution is also computed using another RK method as
    !!
    !!  \f$ Q'_{n+1} = \sum_{i = 1}^{s} b'_{i} \Delta Q_{i} f$
    !!
    !!  An error estimate for the solution advancing scheme is given as \f$ \frac{||Q_{n+ 1} - Q'_{n+1}||}{dt} \f$
    !!
    !!  which is then used to adjust the time step \Delta t by comparing with a user defined tolerance
    !!
    !!  @author Mayank Sharma
    !!  @date   20/1/2017
    !!
    !--------------------------------------------------------------------------------------------------------------
    type, extends(time_integrator_t), public :: embedded_runge_kutta_t


    contains
        procedure   :: iterate

        final       :: destructor
    end type embedded_runge_kutta_t
    !**************************************************************************************************************



    !>  Solve for update dq
    !!
    !!  @author Mayank Sharma
    !!  @date   20/1/2017
    !!
    !--------------------------------------------------------------------------------------------------------------
    subroutine iterate(self,data,nonlinear_solver,linear_solver,preconditioner)
        class(embedded_runge_kutta_t),              intent(inout)   :: self
        type(chidg_data_t),                         intent(inout)   :: data
        class(nonlinear_solver_t),      optional,   intent(inout)   :: nonlinear_solver
        class(linear_solver_t),         optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),        optional,   intent(inout)   :: preconditioner

        !
        ! Variables read from subroutines defining the embedded methods
        !
        integer(ik)                         :: nstage   
        real(rk),            allocatable    :: a(:,:)   
        real(rk),            allocatable    :: b(:,:), b_app(:,:)
        real(rk)                            :: p
        !
        ! Variables defined for use in this module
        !
        real(rk)                            :: t_start, t_end, delta, error_mag
        type(chidgVector_t), allocatable    :: delq(:)
        type(chidgVector_t)                 :: q_n, q_app, q_adv, error   ! q_n, q_app and q_adv are temporary variables
        character(len = :),  allocatable    :: time_scheme
        character(100)                      :: filename
        integer(ik)                         :: itime = 1, nsteps, ielem, wcount, iblk, ieqn, idom, istage, j
        real(rk),            allocatable    :: vals(:)


        wcount = 1
        associate( q => data%sdata%dq, dq => data%sdata%dq, rhs => data%sdata%rhs, lhs => data%sdata%lhs, dt => data%time_manager%dt )

        !
        ! Since embedded methods have a varying time step size, the number of steps isn't known beforehand
        ! The iteration is done through a while loop which quits when t > t_end
        ! TODO: Get from time_manager via chidg.nml?
        !
        ! TODO: WATCH OUT here, time_manager is now part of chidg%data

        t_start = ZERO 
        t_end   = t_start + (dt*self%time_manager%nsteps)
        t       = t_start
    
        !
        ! Get name of the embedded RK method being used
        !
        time_scheme = data%time_manager%get_name()

        !
        ! Get nstage, a, b, and err for a particular method
        ! All methods defined in mod_embedded_RK_methods
        !
        call method_selector(time_scheme,nstage,p,a,b,b_app)

            do 

                !
                ! Copy the solution vector at the present time level to a temporary vector
                !
                q_n = q

                !
                ! Initialize the solution vector used to estimate the error
                ! Initialize the solution vector used to advance the solution
                ! Computed by an accumulated sum in istage loop
                !
                q_app = q_n; q_adv = q_n

                if (allocated(delq)) deallocate(delq)
                allocate(delq(nstage))

                !
                ! Loop through the stages common to both RK methods
                !
                do istage = 1,nstage

                    !
                    ! Update Spatial Residual and Linearization (rhs,lin)
                    !
                    call update_space(data,differentiate=.false.)

                    do idom = 1,data%ndomains()
                        do ielem = 1,data%mesh(idom)%nelem
                            
                            do ieqn = 1,data%eqnset(idom)%prop%nprimary_fields()
                                vals = matmul(data%mesh(idom)%elems(ielem)%invmass, rhs%dom(idom)%vecs(ielem)%getvar(ieqn,itime))
                                call rhs%dom(idom)%vecs(ielem)%setvar(ieqn,itime,vals)
                            end do

                        end do  ! ielem
                    end do  ! idom

                    !
                    ! Calculate stagewise updates and store them in a temporary array
                    !
                    dq = (-dt)*rhs
                    delq(istage) = dq

                    !
                    ! Compute the solution vector update to pass to update_space for the next stage
                    !
                    if (istage < nstage) then
                        do j = 1,istage
                            q = q + (a(istage + 1,j)*delq(j))
                        end do
                    end if

                    !
                    ! Clear Residual and Linearization storage
                    !
                    call lhs%clear()
                    call rhs%clear()

                    !
                    ! Update q_app and q_adv with each stage
                    ! Both methods (order p and p + 1/p - 1) use different update coefficients
                    ! See embedded method definitions in mod_embedded_RK_methods.f90
                    !
                    q_app = add_chidgVector_chidgVector(q_app, (b_app(istage)*delq(istage)))
                    q_adv = add_chidgVector_chidgVector(q_adv, (b(istage)*delq(istage)))

                end do  ! istage

                !
                ! Both q_adv and q_app are chidgVectors
                ! The error estimate is obtained by the L2 vector norm ||q_adv - q_app||
                ! The norm here is calculated across all processors
                !
                error = sub_chidgVector_chidgVector(q_adv - q_app)
                error_mag = error%norm(ChiDG_COMM)/dt
                
                !
                ! Check if the step is acceptable
                ! If the step is not acceptable, dt is adjusted and q is recomputed at the same time level
                ! ttol is a user defined tolerance which is used to judge the error magnitude
                ! TODO: Componentwise tolerance and error using absolute and relative tolerances?
                !

                ! TODO: WATCH OUT here, time_manager is now part of chidg%data and ttol is not a time_manager's object
                
                if (error_mag <= data%time_manager%ttol) then
                    
                    ! Advance time by dt
                    t = t + dt

                    ! Loop exit criteria
                    if (t > t_end) exit

                    ! Accept the solution update
                    q = q_adv



                    ! Print diagnostics
                    call write_line('   R(Q) - Norm:     ', rhs%norm(ChiDG_COMM),delimiter = '')

                    if (wcount == data%time_manager%nwrite) then
                        write(filename,'(I7,A4)') 1000000 + itime, '.plt'
                        call write_tecio_variables_unstructured(data,time(filename),itime + 1)
                        wcount = 0
                    end if

                    itime = itime + 1

                end if
                     
                !
                ! delta is used to adjust the step size
                !
                delta = (HALF*data%time_manager%ttol/error_mag)**(ONE/(p + 1))

                !
                ! Built in safeguards to prevent too big or too small adjustments
                !
                if (delta <= TENTH) then
                    dt = TENTH*dt
                else if (delta >= FOUR)
                    dt = FOUR*dt
                else
                    dt = delta*dt
                end if

                !
                ! Enforce bounds on dt
                ! TODO: dt_max/dt_min should come from input.nml?
                !
                !if (dt >= dt_max) dt = dt_max
                !if (dt <= dt_min) then
                !    print *< 'Minimum dt exceeded'
                !    exit
                !end if

            end do  ! itime

        end associate


    end subroutine iterate
    !**************************************************************************************************************




    !>  
    !!
    !!  @author Mayank Sharma
    !!  @date   20/1/2017
    !!
    !--------------------------------------------------------------------------------------------------------------
    subroutine destructor(self)
        type(embedded_runge_kutta_t),   intent(in)  :: self

    end subroutine destructor
    !**************************************************************************************************************




















end module embedded_runge_kutta
