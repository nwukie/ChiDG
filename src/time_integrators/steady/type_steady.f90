module type_steady
#include <messenger.h>
    use mod_kinds,                      only: rk,ik
    use mod_constants,                  only: ZERO
    use mod_spatial,                    only: update_space
    use mod_chidg_mpi,                  only: GLOBAL_MASTER

    use type_time_integrator_steady,    only: time_integrator_steady_t
    use type_system_assembler,          only: system_assembler_t

    use type_chidg_data,                only: chidg_data_t
    use type_nonlinear_solver,          only: nonlinear_solver_t
    use type_linear_solver,             only: linear_solver_t
    use type_preconditioner,            only: preconditioner_t
    use type_chidg_vector

    use mod_entropy,                    only: compute_entropy_error
    implicit none
    private



    !>  Object implementing solving a system of equations to steady-state.
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/29/2016
    !!
    !----------------------------------------------------------------------------------------
    type, extends(time_integrator_steady_t), public :: steady_t



    contains

        procedure   :: init
        procedure   :: step

    end type steady_t
    !****************************************************************************************





    !>  Object for assembling the implicit system.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2017
    !!
    !----------------------------------------------------------------------------------------
    type, extends(system_assembler_t), public :: assemble_steady_t


    contains

        procedure   :: assemble

    end type assemble_steady_t
    !****************************************************************************************






contains





    !>  Initialize the steady_t time integrator.
    !!
    !!  Main activity here, create the assembler and attach it to the time_integrator 
    !!  object so it can be passed to the nonlinear solver.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2017
    !!
    !----------------------------------------------------------------------------------------
    subroutine init(self)
        class(steady_t),    intent(inout)   :: self

        integer(ik)             :: ierr
        type(assemble_steady_t) :: assemble_steady

        ! Set name
        call self%set_name('Steady')

        ! Allocate assembler
        if (allocated(self%system)) deallocate(self%system)
        allocate(self%system, source=assemble_steady, stat=ierr)
        if (ierr /= 0) call AllocationError

    end subroutine init
    !*****************************************************************************************




    !>  Steady equation set. No time.
    !!
    !!  Solving R(Q) = 0
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/29/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine step(self,data,nonlinear_solver,linear_solver,preconditioner)
        class(steady_t),                        intent(inout)   :: self
        type(chidg_data_t),                     intent(inout)   :: data
        class(nonlinear_solver_t),  optional,   intent(inout)   :: nonlinear_solver
        class(linear_solver_t),     optional,   intent(inout)   :: linear_solver
        class(preconditioner_t),    optional,   intent(inout)   :: preconditioner

        ! Simply solve the nonlinear system. No iteration in time.
        call nonlinear_solver%solve(data,self%system,linear_solver,preconditioner)

        ! Store end residual from nonlinear solver.
        call self%residual_norm%push_back(nonlinear_solver%residual_norm%at(nonlinear_solver%residual_norm%size()))

    end subroutine step
    !******************************************************************************************








    !>  Assemble the system for the 'Steady' equations with no time derivative contributions.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/8/2016
    !!
    !-----------------------------------------------------------------------------------------
    subroutine assemble(self,data,differentiate,timing)
        class(assemble_steady_t),   intent(inout)               :: self
        type(chidg_data_t),         intent(inout)               :: data
        logical,                    intent(in)                  :: differentiate
        real(rk),                   intent(inout),  optional    :: timing

        integer(ik) :: eqn_ID, idom, ielem, itime, imat, ifield, irow_start, icol_start, nterms, nfields, ifield_row, ifield_col
        real(rk),   allocatable     :: elem_field(:), elem_res(:)


        call data%sdata%rhs%clear()
        if (differentiate) call data%sdata%lhs%clear()

        ! Steady equation, so we only need the spatial operators computed.
        data%time_manager%itime = 1
        data%time_manager%t     = ZERO
        call update_space(data,differentiate,timing)


!        !
!        ! Filter q, rhs
!        !
!        do idom = 1,data%mesh%ndomains()
!            eqn_ID = data%mesh%domain(idom)%eqn_ID
!            if (index(trim(data%eqnset(eqn_ID)%name),'Filtered') /= 0) then
!                do ielem = 1,data%mesh%domain(idom)%nelements()
!                    do ifield = 1,data%eqnset(eqn_ID)%prop%nprimary_fields()
!                        
!                        !
!                        ! Get solution field
!                        !
!                        elem_field = data%sdata%q%dom(idom)%vecs(ielem)%getvar(ifield,1)
!                        elem_res   = data%sdata%rhs%dom(idom)%vecs(ielem)%getvar(ifield,1)
!
!                        !
!                        ! Apply spatial filtering
!                        !
!                        if (size(elem_field) > 1) then
!                            if (     index(trim(data%eqnset(eqn_ID)%name),'FP0') /= 0) then
!                                elem_field(2:size(elem_field)) = ZERO
!                                elem_res(  2:size(elem_res))   = ZERO
!                            else if (index(trim(data%eqnset(eqn_ID)%name),'FP1') /= 0) then
!                                elem_field(9:size(elem_field)) = ZERO
!                                elem_res(  9:size(elem_res))   = ZERO
!                            else if (index(trim(data%eqnset(eqn_ID)%name),'FP2') /= 0) then
!                                elem_field(28:size(elem_field)) = ZERO
!                                elem_res(  28:size(elem_res))   = ZERO
!                            else if (index(trim(data%eqnset(eqn_ID)%name),'FP3') /= 0) then
!                                elem_field(65:size(elem_field)) = ZERO
!                                elem_res(  65:size(elem_res))   = ZERO
!                            else if (index(trim(data%eqnset(eqn_ID)%name),'FP4') /= 0) then
!                                elem_field(126:size(elem_field)) = ZERO
!                                elem_res(  126:size(elem_res))   = ZERO
!                            else
!
!                                ! Default, filter to piecewise constant
!                                !   : ZERO all but the constant mode
!                                elem_field(2:size(elem_field)) = ZERO
!                                elem_res(  2:size(elem_res))   = ZERO
!
!                            end if
!
!                            !
!                            ! Store filtered data
!                            !
!                            call data%sdata%q%dom(idom)%vecs(ielem)%setvar(ifield,1,elem_field)
!                            call data%sdata%rhs%dom(idom)%vecs(ielem)%setvar(ifield,1,elem_res)
!
!                        end if
!
!
!
!                    end do !ifield
!                end do !ielem
!
!
!
!
!                !
!                ! Filter lhs: lblks
!                !
!                do ielem = 1,size(data%sdata%lhs%dom(idom)%lblks,1)
!                    do itime = 1,size(data%sdata%lhs%dom(idom)%lblks,2)
!                        do imat = 1,data%sdata%lhs%dom(idom)%lblks(ielem,itime)%size()
!
!                            nterms  = data%sdata%lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%nterms
!                            nfields = data%sdata%lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%nfields
!
!
!                            if (nterms > 1) then
!                                do ifield_row = 1,nfields
!                                    do ifield_col = 1,nfields
!
!                                        irow_start = 1 + nterms*(ifield_row-1)
!                                        icol_start = 1 + nterms*(ifield_col-1)
!
!                                        ! Zero trailing row
!                                        data%sdata%lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(icol_start  ,irow_start+1:irow_start+nterms-1) = ZERO
!                                        ! Zero lower rows
!                                        data%sdata%lhs%dom(idom)%lblks(ielem,itime)%data_(imat)%mat(icol_start+1:icol_start+nterms-1,irow_start:irow_start+nterms-1) = ZERO
!                                    end do !ifield_col
!                                end do !ifield_row
!                            end if
!
!
!                        end do !imat
!                    end do !itime
!                end do !ielem
!        
!
!            end if ! if Filtered
!        end do !idom





        


    end subroutine assemble
    !******************************************************************************************



    
end module type_steady
