module type_linear_solver
    use mod_kinds,          only: rk,ik
    use type_dict,          only: dict_t
    use type_chidg_matrix,  only: chidg_matrix_t
    use type_chidg_vector,  only: chidg_vector_t
    use type_chidg_data,    only: chidg_data_t
    use type_timer,         only: timer_t
    use operator_chidg_mv
    implicit none
    


    !> Abstract type for Matrix Solvers used to implement a common interface
    !! for solving the linear system Ax=b  
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/17/2016
    !!
    !!
    !-----------------------------------------------------------------------------------------------------
    type, public, abstract :: linear_solver_t
        ! Parameters for iterative solver
        real(rk)        :: tol  = 1.e-8_rk              !< Absolute convergence tolerance
        real(rk)        :: rtol = 1.e-10_rk             !< Relative convergence tolerance
        integer(ik)     :: niter = 0                    !< Number of iterations
        integer(ik)     :: nmax  = -1                   !< Max number of iterations, default=-1 (unlimited)
        integer(ik)     :: nkrylov = 2000               !< Size of Krylov subspace
        character(3)    :: orthogonalization = 'CGS'    !< Orthogonalization algorithm: 'CGS', 'MGS'
        integer(ik)     :: silence = 0                  
        real(rk)        :: residual_norm


        ! Parameters for inner iterative solver(preconditioning)
        logical         :: inner_fgmres = .true.
        real(rk)        :: inner_tol  = 1.e-1_rk            !< Absolute convergence tolerance
        real(rk)        :: inner_rtol = 1.e-1_rk            !< Relative convergence tolerance
        integer(ik)     :: inner_nmax = 100                 !< Max number of iterations
        integer(ik)     :: inner_nkrylov = 100              !< Size of Krylov subspace
        integer(ik)     :: inner_silence = -10              !< Silence level of inner iterative method
        character(3)    :: inner_orthogonalization = 'CGS'  !< Orthogonalization algorithm: 'CGS', 'MGS'


        type(timer_t)   :: timer                        !< Timer for linear system solve
        logical         :: report = .true.              !< Flag to enable/disable matrix residual reporting

    contains
    
        procedure   :: init
        procedure   :: tear_down
        procedure   :: set

        procedure(solve_interface), deferred :: solve

        procedure   :: residual
        procedure   :: error

    end type linear_solver_t
    !******************************************************************************************************



    abstract interface
        subroutine solve_interface(self,A,x,b,M,solver_controller,data)
            use type_chidg_matrix,      only: chidg_matrix_t
            use type_chidg_vector,      only: chidg_vector_t
            use type_chidg_data,        only: chidg_data_t
            use type_preconditioner,    only: preconditioner_t
            use type_solver_controller, only: solver_controller_t
            import linear_solver_t

            class(linear_solver_t),     intent(inout)           :: self
            type(chidg_matrix_t),       intent(inout)           :: A
            type(chidg_vector_t),       intent(inout)           :: x
            type(chidg_vector_t),       intent(inout)           :: b
            class(preconditioner_t),    intent(inout), optional :: M
            class(solver_controller_t), intent(inout), optional :: solver_controller
            type(chidg_data_t),         intent(in),    optional :: data
        end subroutine
    end interface



contains


    !> Base initialization for every matrixsolver
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/17/2016
    !!
    !!
    !---------------------------------------------------------------------------------------------------
    subroutine init(self)
        class(linear_solver_t),  intent(inout)   :: self


    end subroutine init
    !***************************************************************************************************






    !> Procedure for setting base matrix_solver options.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/17/2016
    !!
    !!  @param[in]  options     Dictionary containing base matrixsovler options
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine set(self,options)
        class(linear_solver_t),  intent(inout)   :: self
        type(dict_t),           intent(inout)   :: options

        call options%get('tol',               self%tol)
        call options%get('rtol',              self%rtol)
        call options%get('nmax',              self%nmax)
        call options%get('nkrylov',           self%nkrylov)
        call options%get('orthogonalization', self%orthogonalization)

        call options%get('inner_fgmres',            self%inner_fgmres)
        call options%get('inner_tol',               self%inner_tol)
        call options%get('inner_rtol',              self%inner_rtol)
        call options%get('inner_nmax',              self%inner_nmax)
        call options%get('inner_nkrylov',           self%inner_nkrylov)
        call options%get('inner_silence',           self%inner_silence)
        call options%get('inner_orthogonalization', self%inner_orthogonalization)

    end subroutine set
    !***************************************************************************************************




    

    !> Compute the system residual
    !!
    !! Given the system: 
    !!      Ax = b
    !!
    !! The following should be true:
    !!      b - Ax = 0
    !!
    !! So a residual is defined as:
    !!      R = b - Ax
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/17/2016
    !!
    !--------------------------------------------------------------------------------------------------------
    function residual(self,A,x,b) result(r)
        class(linear_solver_t),  intent(in)     :: self
        type(chidg_matrix_t),    intent(inout)  :: A
        type(chidg_vector_t),    intent(inout)  :: x
        type(chidg_vector_t),    intent(in)     :: b

        type(chidg_vector_t) :: r, tmp

        r = x
        call r%clear()

        ! Compute r = b - Ax
        r = b - chidg_mv(A,x)

    end function residual
    !********************************************************************************************************




    !> Compute the L2-norm of the residual vector for the system.
    !!
    !! Given the system: 
    !!      Ax = b
    !!
    !! The following should be true:
    !!      b - Ax = 0
    !!
    !! So a residual is defined as:
    !!      R = b - Ax
    !!
    !! The error metric computed by this function is the L2 norm of the residual:
    !!
    !!  error = ||R||_2
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/17/2016
    !!
    !-------------------------------------------------------------------------------------------------
    function error(self,A,x,b) result(err)
        class(linear_solver_t),  intent(in)     :: self
        type(chidg_matrix_t),    intent(inout)  :: A
        type(chidg_vector_t),    intent(inout)  :: x
        type(chidg_vector_t),    intent(in)     :: b

        type(chidg_vector_t) :: r
        real(rk) :: err

        ! Compute residual
        r = self%residual(A,x,b)

        ! Compute norm
        !err = r%norm()
        err = r%norm(ChiDG_COMM)

    end function error
    !*************************************************************************************************




    !>  Tear down activities.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/25/2019
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine tear_down(self)
        class(linear_solver_t),  intent(inout)  :: self


    end subroutine tear_down
    !*************************************************************************************************

end module type_linear_solver
