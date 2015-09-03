module atype_matrixsolver
    use mod_kinds,      only: rk,ik
    use type_dict,      only: dict_t





    !> Abstract type for Matrix Solvers used to implement a common interface
    !! for solving the linear system Ax=b  
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !-------------------------------------------------------------
    type, public, abstract :: matrixsolver_t

        real(rk)    :: tol = 1.e-8_rk      !< Convergance tolerance for iterative solvers


    contains
    
        procedure   :: init
        procedure   :: set
        procedure(solve_interface), deferred :: solve

    end type matrixsolver_t
    !-------------------------------------------------------------








    abstract interface
        subroutine solve_interface(self,A,x,b)
            use type_blockmatrix,   only: blockmatrix_t
            use type_blockvector,   only: blockvector_t
            import matrixsolver_t

            class(matrixsolver_t),  intent(inout)   :: self
            type(blockmatrix_t),    intent(inout)   :: A
            type(blockvector_t),    intent(inout)   :: x
            type(blockvector_t),    intent(inout)   :: b
        end subroutine
    end interface







contains


    !> Base initialization for every matrixsolver
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !--------------------------------------------------------
    subroutine init(self)
        class(matrixsolver_t),  intent(inout)   :: self


    end subroutine









    !> Procedure for setting base matrix_solver options
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  options     Dictionary containing base matrixsovler options
    !!
    !---------------------------------------------------------
    subroutine set(self,options)
        class(matrixsolver_t),  intent(inout)   :: self
        type(dict_t),           intent(inout)   :: options



        call options%get('tol',self%tol)


    end subroutine
















end module atype_matrixsolver
