module atype_matrixsolver
    use mod_kinds,          only: rk,ik
    use type_dict,          only: dict_t
    use type_blockmatrix,   only: blockmatrix_t
    use type_blockvector
    use type_timer,         only: timer_t
    use operator_mv

    implicit none
    




    !> Abstract type for Matrix Solvers used to implement a common interface
    !! for solving the linear system Ax=b  
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !-------------------------------------------------------------
    type, public, abstract :: matrixsolver_t

        real(rk)        :: tol = 1.e-7_rk      !< Convergance tolerance for iterative solvers

        type(timer_t)   :: timer                !< Timer for linear system solve

        logical         :: report = .true.      !< Flag to enable/disable matrix residual reporting

    contains
    
        procedure   :: init
        procedure   :: set

        procedure(solve_interface), deferred :: solve

        procedure   :: residual
        procedure   :: error

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
    !!
    !!
    !-----------------------------------------------------------
    function residual(self,A,x,b) result(r)
        class(matrixsolver_t),  intent(inout)   :: self
        type(blockmatrix_t),    intent(inout)   :: A
        type(blockvector_t),    intent(inout)   :: x
        type(blockvector_t),    intent(inout)   :: b


        type(blockvector_t) :: r
        real(rk)            :: err
        integer(ik)         :: iparent, ielem, iblk


        r = x
        call r%clear()


        !
        ! Compute Ax
        !
        !do ielem = 1,size(A%lblks,1)
        !    do iblk = 1,size(A%lblks,2)
!
!                if (allocated(A%lblks(ielem,iblk)%mat)) then
!                    iparent = A%lblks(ielem,iblk)%parent()
!                    r%lvecs(ielem)%vec = r%lvecs(ielem)%vec + matmul(A%lblks(ielem,iblk)%mat,x%lvecs(iparent)%vec)
!                end if
!
!            end do ! iblk
!        end do ! ielem



        !
        ! Compute r = b - Ax
        !
        !err_vec = err_vec - b
        r = b - A*x


    end function






















    !> Compute the residual norm
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
    !!
    !-----------------------------------------------------------
    function error(self,A,x,b) result(err)
        class(matrixsolver_t),  intent(inout)   :: self
        type(blockmatrix_t),    intent(inout)   :: A
        type(blockvector_t),    intent(inout)   :: x
        type(blockvector_t),    intent(inout)   :: b


        type(blockvector_t) :: r
        real(rk)            :: err


        !
        ! Allocate residual vector and clear
        !
        !r = x
        !call r%clear()



        !
        ! Compute residual
        !
        r = self%residual(A,x,b)
    

        !
        ! Compute norm
        !
        err = r%norm()


    end function










end module atype_matrixsolver
