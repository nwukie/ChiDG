module type_solver_controller
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: ZERO
    use type_chidg_matrix,      only: chidg_matrix_t
    use type_preconditioner,    only: preconditioner_t
    use mpi_f08,                only: MPI_AllReduce, MPI_LOGICAL, MPI_LOR, mpi_comm
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
    type, public :: solver_controller_t

        ! LHS controller data
        logical :: force_update_lhs            = .false.
        logical :: lhs_updated

        ! Preconditioner controller data
        logical :: force_update_preconditioner = .false.
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
    function update_lhs(self,A,niter,residual_ratio) result(update)
        class(solver_controller_t), intent(inout)   :: self
        type(chidg_matrix_t),       intent(in)      :: A
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
    function update_preconditioner(self,lhs,preconditioner,comm) result(update)
        class(solver_controller_t), intent(inout)   :: self
        type(chidg_matrix_t),       intent(in)      :: lhs
        class(preconditioner_t),    intent(in)      :: preconditioner
        type(mpi_comm),             intent(in)      :: comm

        logical     :: lhs_changed
        logical     :: update, update_global
        integer(ik) :: ierr

        ! Update preconditioner:
        !   1: if lhs has been changed
        !   2: if being forced
        !   3: if preconditioner update was not with this version of the matrix
        update = (any(self%lhs_stamp /= lhs%stamp) .or. &   
                  self%force_update_preconditioner .or. &
                  any(lhs%stamp /= preconditioner%stamp) )

        ! Update last lhs date_and_time
        self%lhs_stamp = lhs%stamp

        ! Turn off forced update
        self%force_update_preconditioner = .false.

        ! Make sure all processes agree on action (some might not have elements, 
        ! so their matrix might not be stamped with an update.)
        call MPI_AllReduce(update,update_global,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
        if (ierr /= 0) call chidg_signal(FATAL,'solver_controller%update_preconditioner: error reducing preconditioner update flat.')

        ! Store global action
        update = update_global

        ! Store action
        self%preconditioner_updated = update

    end function update_preconditioner
    !*************************************************************










end module type_solver_controller
