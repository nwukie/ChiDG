module type_timescheme
    use messenger,          only: write_line
    use mod_kinds,          only: rk,ik
    use atype_matrixsolver, only: matrixsolver_t
    use type_dict,          only: dict_t
    use type_timer,         only: timer_t
    use type_rvector,       only: rvector_t
    use type_ivector,       only: ivector_t
    use type_chidg_data,    only: chidg_data_t
    implicit none


    !> solver abstract type definition
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !-------------------------------------------------------------------------------
    type, abstract, public  :: timescheme_t

        real(rk)        :: testing
        logical         :: solverInitialized = .false.



        ! OPTIONS
        real(rk)        :: cfl0     = 1.0_rk        !< Initial CFL number
        real(rk)        :: dt       = 0.001_rk      !< Time-step increment
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

        ! Must define these procedures in the extended type
        procedure(data_interface),   deferred   :: solve

        procedure   :: set
        procedure   :: report
    end type timescheme_t
    !------------------------------------------------------------------------------









    !==================================================
    !
    !   solver deferred procedure interfaces
    !
    !==================================================
    
    abstract interface
        subroutine self_interface(self)
            import timescheme_t
            class(timescheme_t),   intent(inout)   :: self
        end subroutine
    end interface



    abstract interface
        subroutine init_interface(self,data,options)
            use type_chidg_data, only: chidg_data_t
            use type_dict,       only: dict_t
            import timescheme_t
            class(timescheme_t),    intent(inout)   :: self
            type(chidg_data_t),     intent(inout)   :: data
            type(dict_t), optional, intent(inout)   :: options
        end subroutine
    end interface





    ! Interface for passing a domain_t type
    abstract interface
        subroutine data_interface(self,data,matrixsolver,preconditioner)
            use type_chidg_data,        only: chidg_data_t
            use atype_matrixsolver,     only: matrixsolver_t
            use type_preconditioner,    only: preconditioner_t
            import timescheme_t
            class(timescheme_t),                    intent(inout)   :: self
            type(chidg_data_t),                     intent(inout)   :: data
            class(matrixsolver_t),      optional,   intent(inout)   :: matrixsolver
            class(preconditioner_t),    optional,   intent(inout)   :: preconditioner
        end subroutine
    end interface




contains

    !> Common time_scheme initialization interface.
    !!      - Call initialization for options if present
    !!      - Call user-specified initialization routine for concrete type
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[inout]   domains     Array of domains
    !!  @param[inout]   options     Dictionary containing options
    !-------------------------------------------------------------------------
    subroutine init(self,data)
        class(timescheme_t),    intent(inout)   :: self
        type(chidg_data_t),     intent(inout)   :: data


        !
        ! Call user-specified initialization
        !
        call self%init_spec(data)



        self%solverInitialized = .true.
    end subroutine






    !> Procedure for setting base time_scheme options
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  options     Dictionary containing base solver options
    !!
    !------------------------------------------------------------------------
    subroutine set(self,options)
        class(timescheme_t),    intent(inout)   :: self
        type(dict_t),           intent(inout)   :: options


        call options%get('dt',self%dt)
        call options%get('tol',self%tol)
        call options%get('nsteps',self%nsteps)
        call options%get('nwrite',self%nwrite)
        call options%get('cfl0',self%cfl0)
    end subroutine







    !> Default blank initialization-specialization routine.
    !! This can be overwritten with specific instructions for a conrete
    !! time_scheme.
    !!
    !!  @author Nathan A. Wukie
    !!
    !-----------------------------------------------------------------------
    subroutine init_spec(self,data,options)
        class(timescheme_t),    intent(inout)   :: self
        type(chidg_data_t),     intent(inout)   :: data
        type(dict_t), optional, intent(inout)   :: options



    end subroutine










    !>
    !!
    !!
    !!
    !!
    !!
    !!
    !---------------------------------------------------------------------------------------------------------------------
    subroutine report(self)
        class(timescheme_t),   intent(in)  :: self

        integer(ik) :: i

        real(rk)    :: residual_time, residual_norm, matrix_time, total_residual, total_matrix
        integer(ik) :: matrix_iterations


        !
        ! Time scheme header
        !
        call write_line(' ')
        call write_line('---------------------------------   Time Scheme Report  ----------------------------------')
        call write_line('Newton iterations: ', self%newton_iterations%at(1))
        call write_line('Total time: ', self%total_time%at(1))

        call write_line(' ')
        call write_line('------------------------------------------------------------------------------------------')



        !
        ! Print per-iteration report
        !
        call write_line('        Residual time', '              Norm[R]', '             Matrix time', '      Matrix iterations')

        ! Loop through stored data and print for each newton iteration
        do i = 1,self%residual_time%size()
            residual_time     = self%residual_time%at(i)
            residual_norm     = self%residual_norm%at(i)
            matrix_time       = self%matrix_time%at(i)
            matrix_iterations = self%matrix_iterations%at(i)
            
            call write_line(residual_time, residual_norm, matrix_time, matrix_iterations, delimiter=', ')
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
        call write_line('Total residual time: ', total_residual)
        call write_line('Total matrix time  : ', total_matrix)



        call write_line('------------------------------------------------------------------------------------------')


    end subroutine report
    !###############################################################################################################








end module type_timescheme
