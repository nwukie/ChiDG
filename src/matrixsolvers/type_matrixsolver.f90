module type_matrixsolver
    use mod_kinds,          only: rk,ik
    use type_dict,          only: dict_t
    use type_chidgMatrix,   only: chidgMatrix_t
    use type_chidgVector
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
    type, public, abstract :: matrixsolver_t

!        character(len=:), allocatable   :: name_

        ! OPTIONS
        real(rk)        :: tol   = 1.e-8_rk     !< Convergance tolerance for iterative solvers
        integer(ik)     :: niter = 0


        type(timer_t)   :: timer                !< Timer for linear system solve
        logical         :: report = .true.      !< Flag to enable/disable matrix residual reporting

    contains
    
        procedure   :: init
!        procedure   :: set_options
        procedure   :: set
!        procedure   :: get_name

        procedure(solve_interface), deferred :: solve

        procedure   :: residual
        procedure   :: error

    end type matrixsolver_t
    !******************************************************************************************************








    abstract interface
        subroutine solve_interface(self,A,x,b,M)
            use type_chidgMatrix,       only: chidgMatrix_t
            use type_chidgVector,       only: chidgVector_t
            use type_preconditioner,    only: preconditioner_t
            import matrixsolver_t

            class(matrixsolver_t),      intent(inout)           :: self
            type(chidgMatrix_t),        intent(inout)           :: A
            type(chidgVector_t),        intent(inout)           :: x
            type(chidgVector_t),        intent(inout)           :: b
            class(preconditioner_t),    intent(inout), optional :: M
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
        class(matrixsolver_t),  intent(inout)   :: self


    end subroutine init
    !***************************************************************************************************







!    !>
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   2/17/2016
!    !!
!    !!
!    !!
!    !------------------------------------------------------------------------------------------
!    subroutine set_options(self)
!        class(matrixsolver_t),  intent(inout)   :: self
!
!
!
!    end subroutine set_options
!    !*******************************************************************************************









    !> Procedure for setting base matrix_solver options.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/17/2016
    !!
    !!  @param[in]  options     Dictionary containing base matrixsovler options
    !!
    !--------------------------------------------------------------------------------------------------
    subroutine set(self,options)
        class(matrixsolver_t),  intent(inout)   :: self
        type(dict_t),           intent(inout)   :: options

        call options%get('tol',     self%tol)
        !call options%get('niter',   self%niter)

    end subroutine set
    !***************************************************************************************************





!    !>  Return the name of the matrix solver.
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   2/17/2016
!    !!
!    !!
!    !!
!    !-----------------------------------------------------------------------------------------
!    function get_name(self) result(res)
!        class(matrixsolver_t),  intent(in)  :: self
!
!        character(len=:), allocatable :: res
!
!        res = trim(self%name_)
!
!    end function get_name
!    !******************************************************************************************







    

    !> Compute the system residual
    !!
    !! Given the system: 
    !!      Ax = b
    !!
    !!
    !! The following should be true:
    !!      b - Ax = 0
    !!
    !!
    !! So a residual is defined as:
    !!      R = b - Ax
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/17/2016
    !!
    !--------------------------------------------------------------------------------------------------------
    function residual(self,A,x,b) result(r)
        class(matrixsolver_t),  intent(inout)   :: self
        type(chidgMatrix_t),    intent(inout)   :: A
        type(chidgVector_t),    intent(inout)   :: x
        type(chidgVector_t),    intent(inout)   :: b


        type(chidgVector_t) :: r
        real(rk)            :: err
        integer(ik)         :: iparent, ielem, iblk


        r = x
        call r%clear()


        !
        ! Compute r = b - Ax
        !
        r = b - A*x


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
        class(matrixsolver_t),  intent(inout)   :: self
        type(chidgMatrix_t),    intent(inout)   :: A
        type(chidgVector_t),    intent(inout)   :: x
        type(chidgVector_t),    intent(inout)   :: b


        type(chidgVector_t) :: r
        real(rk)            :: err



        !
        ! Compute residual
        !
        r = self%residual(A,x,b)
    

        !
        ! Compute norm
        !
        err = r%norm()


    end function error
    !*************************************************************************************************










end module type_matrixsolver
