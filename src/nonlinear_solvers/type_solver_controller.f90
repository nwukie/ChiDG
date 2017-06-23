module type_solver_controller
    use mod_kinds,          only: rk, ik
    use mod_constants,      only: ZERO
    use type_chidg_matrix,  only: chidg_matrix_t
    implicit none


    !>  A class for controlling the lhs update inside of 
    !!  a nonlinear solver.
    !!
    !!
    !!  Default behaviors:
    !!  ------------------
    !!  update_lhs:             Always
    !!  update_preconditioner:  If matrix has changed
    !!
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/21/2017
    !!
    !---------------------------------------------------------------
    type, abstract, public :: solver_controller_t

        ! LHS controller data
        logical :: force_update_lhs            = .true.
        logical :: lhs_updated

        ! Preconditioner controller data
        logical :: force_update_preconditioner = .true.
        integer :: lhs_stamp(8) = ZERO                      ! lhs date_time_stamp to recognize when lhs has been updated
        logical :: preconditioner_updated      

    contains

        procedure   :: update_lhs
        procedure   :: update_preconditioner

    end type solver_controller_t
    !***************************************************************


contains




    !>  Control update of the lhs matrix.
    !!
    !!  Default behavior: always update lhs
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/22/2017
    !!
    !-------------------------------------------------------------
    function update_lhs(self,niter,residual_ratio) result(update)
        class(solver_controller_t), intent(inout)   :: self
        integer(ik),                intent(in)      :: niter
        real(rk),                   intent(in)      :: residual_ratio

        logical :: update 
        
        ! Default: update lhs
        update = .true.

        ! Store action
        self%lhs_updated = update

    end function update_lhs
    !*************************************************************





    !>  Control update of the preconditioner.
    !!
    !!  Default behavior: update preconditioner if matrix has changed.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   6/22/2017
    !!
    !-------------------------------------------------------------
    function update_preconditioner(self,lhs) result(update)
        class(solver_controller_t), intent(inout)   :: self
        type(chidg_matrix_t),       intent(inout)   :: lhs

        logical :: lhs_changed
        logical :: update

        ! Update preconditioner:
        !   1: if lhs has been changed
        !   2: if being forced
        update = (any(self%lhs_stamp /= lhs%stamp) .or. self%force_update_preconditioner)
        
        ! Update last lhs date_and_time
        self%lhs_stamp = lhs%stamp

        ! Turn off forced update
        self%force_update_preconditioner = .false.

        ! Store action
        self%preconditioner_updated = update

    end function update_preconditioner
    !*************************************************************










end module type_solver_controller
