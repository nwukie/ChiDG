module type_nonlinear_solver
    use messenger,          only: write_line
    use mod_kinds,          only: rk,ik
    use type_linear_solver, only: linear_solver_t
    use type_dict,          only: dict_t
    use type_timer,         only: timer_t
    use type_rvector,       only: rvector_t
    use type_ivector,       only: ivector_t
    use type_chidg_data,    only: chidg_data_t
    implicit none


    !> solver abstract type definition
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !------------------------------------------------------------------------------------------------------
    type, abstract, public  :: nonlinear_solver_t

        logical         :: solverInitialized = .false.

        ! OPTIONS
        real(rk)        :: cfl0     = 1.0_rk        !< Initial CFL number
        real(rk)        :: tol      = 1.e-13_rk     !< Convergence tolerance
        integer(ik)     :: nsteps   = 100           !< Number of time steps to compute
        integer(ik)     :: nwrite   = 10            !< Write data every 'nwrite' steps


        type(timer_t)   :: timer                    !< Timer data-type


        ! Data logs
        type(rvector_t) :: residual_norm
        type(rvector_t) :: residual_time
        type(ivector_t) :: matrix_iterations
        type(rvector_t) :: matrix_time
        type(ivector_t) :: newton_iterations
        type(rvector_t) :: total_time


    contains

        procedure   :: init
        procedure   :: init_spec

        procedure(data_interface),   deferred   :: solve    ! Must define this procedures in the extended type

        procedure   :: set
        procedure   :: report

    end type nonlinear_solver_t
    !******************************************************************************************************









    !==================================================
    !
    !   solver deferred procedure interfaces
    !
    !==================================================
    
    abstract interface
        subroutine self_interface(self)
            import nonlinear_solver_t
            class(nonlinear_solver_t),   intent(inout)   :: self
        end subroutine
    end interface



    abstract interface
        subroutine init_interface(self,data,options)
            use type_chidg_data, only: chidg_data_t
            use type_dict,       only: dict_t
            import nonlinear_solver_t
            class(nonlinear_solver_t),  intent(inout)   :: self
            type(chidg_data_t),         intent(inout)   :: data
            type(dict_t), optional,     intent(inout)   :: options
        end subroutine
    end interface





    ! Interface for passing a domain_t type
    abstract interface
        subroutine data_interface(self,data,linear_solver,preconditioner)
            use type_chidg_data,        only: chidg_data_t
            use type_linear_solver,     only: linear_solver_t
            use type_preconditioner,    only: preconditioner_t
            import nonlinear_solver_t
            class(nonlinear_solver_t),              intent(inout)   :: self
            type(chidg_data_t),                     intent(inout)   :: data
            class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
            class(preconditioner_t),    optional,   intent(inout)   :: preconditioner
        end subroutine
    end interface




contains




    !>  Common time_scheme initialization interface.
    !!      - Call initialization for options if present
    !!      - Call user-specified initialization routine for concrete type
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!  @param[inout]   domains     Array of domains
    !!  @param[inout]   options     Dictionary containing options
    !!
    !-------------------------------------------------------------------------------------------------------------
    subroutine init(self,data)
        class(nonlinear_solver_t),  intent(inout)   :: self
        type(chidg_data_t),         intent(inout)   :: data


        !
        ! Call user-specified initialization
        !
        call self%init_spec(data)


        self%solverInitialized = .true.

    end subroutine init
    !************************************************************************************************************










    !> Procedure for setting base time_scheme options
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!  @param[in]  options     Dictionary containing base solver options
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine set(self,options)
        class(nonlinear_solver_t),  intent(inout)   :: self
        type(dict_t),               intent(inout)   :: options


!        call options%get('tol',self%tol)
!        call options%get('nsteps',self%nsteps)
!        call options%get('nwrite',self%nwrite)

    end subroutine set
    !*************************************************************************************************************












    !> Default blank initialization-specialization routine.
    !! This can be overwritten with specific instructions for a conrete
    !! time_scheme.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !-------------------------------------------------------------------------------------------------------------
    subroutine init_spec(self,data,options)
        class(nonlinear_solver_t),  intent(inout)   :: self
        type(chidg_data_t),         intent(inout)   :: data
        type(dict_t), optional,     intent(inout)   :: options



    end subroutine init_spec
    !*************************************************************************************************************










    !>  Print timescheme report
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------------------------------
    subroutine report(self)
        class(nonlinear_solver_t),   intent(in)  :: self

        integer(ik) :: i

        real(rk)    :: residual_time, residual_norm, matrix_time, total_residual, total_matrix
        integer(ik) :: matrix_iterations


        !
        ! Time scheme header
        !
        call write_line(' ')
        call write_line('---------------------------------   Time Scheme Report  ----------------------------------')
        call write_line('Newton iterations: ', self%newton_iterations%at(1), columns=.True., column_width=20)
        call write_line('Total time: ', self%total_time%at(1), columns=.True., column_width=20)
        call write_line(' ')
        call write_line('------------------------------------------------------------------------------------------')



        !
        ! Print per-iteration report
        !
        call write_line('Residual time', 'Norm[R]', 'Matrix time', 'Matrix iterations', columns=.True., column_width=20)


        !
        ! Loop through stored data and print for each newton iteration
        !
        do i = 1,self%residual_time%size()
            residual_time     = self%residual_time%at(i)
            residual_norm     = self%residual_norm%at(i)
            matrix_time       = self%matrix_time%at(i)
            matrix_iterations = self%matrix_iterations%at(i)
            
            call write_line(residual_time, residual_norm, matrix_time, matrix_iterations, delimiter=', ', columns=.True., column_width=20)
        end do


        !
        ! Accumulate total residual and matrix solver compute times
        !
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
    !***************************************************************************************************************








end module type_nonlinear_solver
