module type_time_integrator
#include <messenger.h>
    use messenger,              only: write_line
    use mod_kinds,              only: rk,ik
    use type_preconditioner,    only: preconditioner_t
    use type_linear_solver,     only: linear_solver_t
    use type_nonlinear_solver,  only: nonlinear_solver_t
    use type_dict,              only: dict_t
    use type_timer,             only: timer_t
    use type_rvector,           only: rvector_t
    use type_ivector,           only: ivector_t
    use type_chidg_data,        only: chidg_data_t
    use type_system_assembler,  only: system_assembler_t
    implicit none


    !>  Abstraction for time integrators.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!  @date   2/7/2017
    !!
    !---------------------------------------------------------------------------------------
    type, abstract, public  :: time_integrator_t

        character(:),   allocatable :: name_

        class(system_assembler_t),  allocatable :: system


        ! Data logs
        type(timer_t)       :: timer
        type(rvector_t)     :: residual_norm
        type(rvector_t)     :: residual_time
        type(ivector_t)     :: matrix_iterations
        type(rvector_t)     :: matrix_time
        type(ivector_t)     :: newton_iterations
        type(rvector_t)     :: total_time

        logical             :: solverInitialized = .false.

    contains

        procedure   :: init
        procedure   :: set
        procedure   :: report
        procedure   :: set_name
        procedure   :: get_name

        ! Adjoint infrastructure
        procedure   :: compute_functionals
        procedure   :: compute_adjoint

        ! Must define this procedure in any extended type
        procedure(step_interface),   deferred   :: step                 ! define the time integrator stepping procedure
        procedure(state_interface),  deferred   :: initialize_state     ! process a read solution for the current time-integrator
        procedure(write_interface),  deferred   :: write_time_options   ! write time integrator options to .h5 file
        procedure(read_interface),   deferred   :: read_time_options    ! read time integrator options from .h5 file
        procedure(post_interface),   deferred   :: process_data_for_output   

    end type time_integrator_t
    !***************************************************************************************









    !==================================================
    !
    !   solver deferred procedure interfaces
    !
    !==================================================
    
    ! Interface for passing a domain_t type
    abstract interface
        subroutine step_interface(self,data,nonlinear_solver,linear_solver,preconditioner)
            use type_chidg_data,        only: chidg_data_t
            use type_nonlinear_solver,  only: nonlinear_solver_t
            use type_linear_solver,     only: linear_solver_t
            use type_preconditioner,    only: preconditioner_t
            import time_integrator_t
            class(time_integrator_t),               intent(inout)   :: self
            type(chidg_data_t),                     intent(inout)   :: data
            class(nonlinear_solver_t),  optional,   intent(inout)   :: nonlinear_solver
            class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
            class(preconditioner_t),    optional,   intent(inout)   :: preconditioner
        end subroutine
    end interface


    ! Interface for passing a domain_t type
    abstract interface
        subroutine state_interface(self,data)
            use type_chidg_data,        only: chidg_data_t
            import time_integrator_t
            class(time_integrator_t),   intent(inout)   :: self
            type(chidg_data_t),         intent(inout)   :: data
        end subroutine
    end interface


    ! Interface for writing time_integrator options to hdf5 file
    abstract interface
        subroutine write_interface(self,data,filename)
            use type_chidg_data,        only: chidg_data_t
            import time_integrator_t
            class(time_integrator_t),   intent(inout)   :: self
            type(chidg_data_t),         intent(inout)   :: data
            character(*),               intent(in)      :: filename
        end subroutine
    end interface


    ! Interface for reading time_integrator options from .h5 file
    abstract interface
        subroutine read_interface(self,data,filename,read_type)
            use type_chidg_data,        only: chidg_data_t
            import time_integrator_t
            class(time_integrator_t),   intent(inout)   :: self
            type(chidg_data_t),         intent(inout)   :: data
            character(*),               intent(in)      :: filename
            character(*),               intent(in)      :: read_type
        end subroutine
    end interface


    ! Interface for modifying data for post processing
    abstract interface
        subroutine post_interface(self,data)
            use type_chidg_data,        only: chidg_data_t
            import time_integrator_t
            class(time_integrator_t),   intent(inout)   :: self
            type(chidg_data_t),         intent(inout)   :: data
        end subroutine
    end interface

!    ! Interface for modifying data for post processing
!    abstract interface
!        subroutine adj_interface(self,data,linear_solver,preconditioner)
!            use type_chidg_data,        only: chidg_data_t
!            use type_linear_solver,     only: linear_solver_t
!            use type_preconditioner,    only: preconditioner_t
!            import time_integrator_t
!            class(time_integrator_t),               intent(inout)   :: self
!            type(chidg_data_t),                     intent(inout)   :: data
!            class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
!            class(preconditioner_t),    optional,   intent(inout)   :: preconditioner
!        end subroutine
!    end interface


contains




    !>  Common time_integrator initialization interface.
    !!      - Call initialization for options if present
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !---------------------------------------------------------------------------------
    subroutine init(self)
        class(time_integrator_t),   intent(inout)   :: self

        self%solverInitialized = .true.

    end subroutine init
    !*********************************************************************************


    !>  Set name
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/4/2017
    !!
    !---------------------------------------------------------------------------------
    subroutine set_name(self,name_)
        class(time_integrator_t),   intent(inout)   :: self
        character(*),               intent(in)      :: name_

        self%name_ = trim(name_)

    end subroutine set_name
    !*********************************************************************************



    !>  Return self%name_
    !!
    !!  @author Nathan A. Wukie
    !!  @date   12/4/2017
    !!
    !---------------------------------------------------------------------------------
    function get_name(self) result(name_)
        class(time_integrator_t),   intent(in)  :: self

        character(:),   allocatable :: name_

        name_ = trim(self%name_)

    end function get_name
    !*********************************************************************************



    !> Procedure for setting base time_integrator options
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!  @param[in]  options     Dictionary containing base solver options
    !!
    !! NOTE: time options (toptions) are not used anymore and time info are stored 
    !!       by time_manager
    !!
    !--------------------------------------------------------------------------------
    subroutine set(self,options)
        class(time_integrator_t),   intent(inout)   :: self
        type(dict_t),               intent(inout)   :: options


        !call options%get('dt',self%dt)
        !call options%get('nsteps',self%nsteps)
        !call options%get('nwrite',self%nwrite)
        !call options%get('tol',self%tol)
        !call options%get('cfl0',self%cfl0)

    end subroutine set
    !******************************************************************************************







    !> Default blank initialization-specialization routine.
    !! This can be overwritten with specific instructions for a conrete
    !! time_integrator.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !------------------------------------------------------------------------------------------
    subroutine init_spec(self,data,options)
        class(time_integrator_t),   intent(inout)   :: self
        type(chidg_data_t),         intent(inout)   :: data
        type(dict_t),   optional,   intent(inout)   :: options



    end subroutine init_spec
    !******************************************************************************************



    !>  Computation of functionals
    !!
    !!  @author Matteo Ugolotti
    !!  @date   6/17/2017
    !!
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute_functionals(self,data)
        class(time_integrator_t),               intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        



    end subroutine compute_functionals
    !******************************************************************************************




    !>  Compute adjoint.
    !!
    !!  @author Matteo Ugolotti
    !!  @date   6/17/2017
    !!
    !-----------------------------------------------------------------------------------------
    subroutine compute_adjoint(self,data,linear_solver,preconditioner)
        class(time_integrator_t),               intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner
        
        call chidg_signal_one(FATAL,"compute_adjoint: adjoint method is not yet implemented for this time integrator ",data%time_manager%get_name())

    end subroutine compute_adjoint
    !******************************************************************************************



    !>  Print timeintegrator report
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine report(self)
        class(time_integrator_t),   intent(in)  :: self

        integer(ik) :: i

        real(rk)    :: residual_time, residual_norm, matrix_time, total_residual, total_matrix
        integer(ik) :: matrix_iterations


        ! Time integrator header
        call write_line(' ')
        call write_line('---------------------------------   Time Integrator Report  ----------------------------------')
        call write_line('Newton iterations: ', self%newton_iterations%at(1), columns=.True., column_width=20)
        call write_line('Total time: ', self%total_time%at(1), columns=.True., column_width=20)
        call write_line(' ')
        call write_line('------------------------------------------------------------------------------------------')


        ! Print per-iteration report
        call write_line('Residual time', 'Norm[R]', 'Matrix time', 'Matrix iterations', columns=.True., column_width=20)

        ! Loop through stored data and print for each newton iteration
        do i = 1,self%residual_time%size()
            residual_time     = self%residual_time%at(i)
            residual_norm     = self%residual_norm%at(i)
            matrix_time       = self%matrix_time%at(i)
            matrix_iterations = self%matrix_iterations%at(i)
            
            call write_line(residual_time, residual_norm, matrix_time, matrix_iterations, delimiter=', ', columns=.True., column_width=20)
        end do


        ! Accumulate total residual and matrix solver compute times
        total_residual = 0._rk
        total_matrix   = 0._rk
        do i = 1,self%residual_time%size()
            total_residual = total_residual + self%residual_time%at(i)
            total_matrix   = total_matrix   + self%matrix_time%at(i)
        end do



        call write_line(' ')
        call write_line('Total residual time: ', total_residual, columns=.True., column_width=20)
        call write_line('Total matrix time  : ', total_matrix,   columns=.True., column_width=20)
        call write_line('------------------------------------------------------------------------------------------')


    end subroutine report
    !******************************************************************************************








end module type_time_integrator
