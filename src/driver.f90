!>  Chimera-based, discontinuous Galerkin equation solver
!!
!!  This program is designed to solve partial differential equations,
!!  and systems of partial differential equations, using the discontinuous
!!  Galerkin method for spatial discretization using Chimera, overset grids to
!!  represent the simulation domain.
!!
!!  @author Nathan A. Wukie
!!  @date   1/31/2016
!!
!!
!---------------------------------------------------------------------------------------------
program driver
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use type_chidg,                 only: chidg_t
    use type_chidg_manager,         only: chidg_manager_t
    use type_function,              only: function_t
    use mod_function,               only: create_function
    use mod_chidg_mpi,              only: GLOBAL_MASTER, ChiDG_COMM
    use mod_io

    ! Actions
    use mod_chidg_edit,         only: chidg_edit
    use mod_chidg_convert,      only: chidg_convert
    use mod_chidg_post,         only: chidg_post,chidg_post_vtk

    
    !
    ! Variable declarations
    !
    implicit none
    type(chidg_manager_t)                       :: manager
    type(chidg_t)                               :: chidg


    integer                                     :: narg, iorder, ierr
    character(len=1024)                         :: chidg_action, filename
    class(function_t),              allocatable :: constant, monopole, fcn, polynomial





    !
    ! Check for command-line arguments
    !
    narg = command_argument_count()


    !
    ! Execute ChiDG calculation
    !
    if ( narg == 0 ) then



        !
        ! Initialize ChiDG environment
        !
        call chidg%start_up('mpi')
        call chidg%start_up('namelist')
        call chidg%start_up('core')



        !
        ! Set ChiDG Algorithms
        !
!        call chidg%set('Time Integrator' , algorithm=time_integrator,  options=toptions)
        call chidg%set('Time Integrator' , algorithm=time_integrator)
        call chidg%set('Nonlinear Solver', algorithm=nonlinear_solver, options=noptions)
        call chidg%set('Linear Solver'   , algorithm=linear_solver,    options=loptions)
        call chidg%set('Preconditioner'  , algorithm=preconditioner                    )


        !
        ! Set ChiDG Files, Order, etc.
        !
        !call chidg%set('Grid File',         file=grid_file              )
        !call chidg%set('Solution File In',  file=solution_file_in       )
        !call chidg%set('Solution File Out', file=solution_file_out      )
        call chidg%set('Solution Order',    integer_input=solution_order)



        !
        ! Read grid and boundary condition data
        !
        call chidg%read_grid(gridfile)
        call chidg%read_boundaryconditions(gridfile)



        !
        ! Initialize communication, storage, auxiliary fields
        !
        call manager%process(chidg)
        call chidg%init('all')



        !
        ! Initialize solution
        !
        if (solutionfile_in == 'none') then

            !
            ! Set initial solution
            !
!            call create_function(fcn,'gaussian')
!            call fcn%set_option('b_x',0._rk)
!            call fcn%set_option('b_y',3.14_rk)
!            call fcn%set_option('b_z',0._rk)
!            call fcn%set_option('c',0.5_rk)
!            call chidg%data%sdata%q_in%project(chidg%data%mesh,fcn,1)
!            call create_function(constant,'constant')
!            call constant%set_option('val',0._rk)
!            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,1)


!            call polynomial%set_option('f',3.5_rk)
!            call create_function(polynomial,'polynomial')
!
!            ! d
!            call create_function(constant,'constant')
!            call constant%set_option('val',0.001_rk)
!            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,1)


            call create_function(constant,'constant')

            ! rho
            call constant%set_option('val',1.15_rk)
            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,1)

            ! rho_u
            call constant%set_option('val',0.0_rk)
            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,2)

            ! rho_v
            call constant%set_option('val',0.0_rk)
            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,3)

            ! rho_w
            call constant%set_option('val',50.0_rk)
            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,4)

            ! rho_E
            call constant%set_option('val',248000.0_rk)
            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,5)

!            ! rho_nutilde
!            call constant%set_option('val',0.00003_rk)
!            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,6)
!
!            ! eps
!            call constant%set_option('val',0.000001_rk)
!            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,7)

        else

            !
            ! TODO: put in check that solutionfile actually contains solution
            !
            call chidg%read_solution(solutionfile_in)

        end if



        !
        ! Run ChiDG simulation
        !
        call chidg%report('before')

        call chidg%run(write_initial=initial_write, write_final=final_write)

        call chidg%report('after')





        !
        ! Close ChiDG
        !
        call chidg%shut_down('core')
        call chidg%shut_down('mpi')









    !
    ! ChiDG tool execution. 2 arguments.
    !
    else if ( narg == 2 ) then


        call get_command_argument(1,chidg_action)
        call get_command_argument(2,filename)
        chidg_action = trim(chidg_action)
        filename = trim(filename)
        

        !
        ! Initialize ChiDG environment
        !
        call chidg%start_up('core')


        !
        ! Select ChiDG action
        !
        if ( trim(chidg_action) == 'edit' ) then
            call chidg_edit(trim(filename))

        else if ( trim(chidg_action) == 'convert' ) then
            call chidg_convert(trim(filename))

        else if ( trim(chidg_action) == 'post' ) then
            call chidg_post(trim(filename))
            call chidg_post_vtk(trim(filename))

        else
            call chidg_signal(FATAL,"chidg: unrecognized action '"//trim(chidg_action)//"'. Valid options are: 'edit', 'convert'")

        end if


        !
        ! Close ChiDG interface
        !
        call chidg%shut_down('core')


    else
        call chidg_signal(FATAL,"chidg: invalid number of arguments. Expecting (0) arguments: 'chidg'. or (2) arguments: 'chidg action file'.")
    end if












end program driver
