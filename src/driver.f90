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
    use mod_chidg_mpi,              only: GLOBAL_MASTER, ChiDG_COMM, IRANK
    use mod_io

    ! Actions
    use mod_chidg_edit,         only: chidg_edit
    use mod_chidg_convert,      only: chidg_convert
    use mod_chidg_post,         only: chidg_post, chidg_post_vtk, chidg_post_matplotlib
    use mod_chidg_airfoil,      only: chidg_airfoil
    use mod_chidg_clone,        only: chidg_clone

    
    !
    ! Variable declarations
    !
    implicit none
    type(chidg_manager_t)                       :: manager
    type(chidg_t)                               :: chidg


    integer                                     :: narg, iorder, ierr
    character(len=1024)                         :: chidg_action, filename, grid_file, solution_file, file_a, file_b
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
        call chidg%read_mesh(gridfile)

        !
        ! Initialize communication, storage, auxiliary fields
        !
        call manager%process(chidg)



        !
        ! Initialize solution
        !
        if (solutionfile_in == 'none') then

            !
            ! Set initial solution
            !
!            call create_function(fcn,'gaussian')
!            call fcn%set_option('b_x',1._rk)
!            call fcn%set_option('b_y',1._rk)
!            call fcn%set_option('b_z',1._rk)
!            call fcn%set_option('c',0.5_rk)
!            call chidg%data%sdata%q_in%project(chidg%data%mesh,fcn,1)
!            call create_function(constant,'constant')
!            call constant%set_option('val',0._rk)
!            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,1)


!            call polynomial%set_option('f',3.5_rk)
!            call create_function(polynomial,'polynomial')


!            ! d
!            !call create_function(constant,'constant')
!            !call constant%set_option('val',10000.0_rk)
!            call create_function(constant,'Radius')
!            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,1)

            


            call create_function(constant,'constant')


            ! rho
            call constant%set_option('val',1.1_rk)
            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,1)

            ! rho_u
            call constant%set_option('val',10.0_rk)
            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,2)

            ! rho_v
            call constant%set_option('val',0.0_rk)
            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,3)

            ! rho_w
            call constant%set_option('val',0.0_rk)
            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,4)

            ! rho_E
            call constant%set_option('val',240000._rk)
            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,5)

!            ! rho_nutilde
!            call constant%set_option('val',0.00009_rk)
!            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,6)
!
!            ! eps
!            call constant%set_option('val',0.000001_rk)
!            call chidg%data%sdata%q_in%project(chidg%data%mesh,constant,7)

        else

            !
            ! TODO: put in check that solutionfile actually contains solution
            !
            call chidg%read_fields(solutionfile_in)

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
    ! Check if executing 'action'
    !
    else if ( narg > 1 ) then

        ! Get 'action'
        call get_command_argument(1,chidg_action)
        call chidg%start_up('core')

        !
        ! Select 'action'
        ! 
        select case (trim(chidg_action))
            case ('edit')
                if (narg /= 2) call chidg_signal(FATAL,"The 'edit' action expects: chidg edit filename.h5")
                call get_command_argument(2,filename)
                call chidg_edit(trim(filename))

            case ('convert')
                if (narg /= 2) call chidg_signal(FATAL,"The 'convert' action expects: chidg convert filename.x")
                call get_command_argument(2,filename)
                call chidg_convert(trim(filename))

            case ('post')
                if (narg /= 3) call chidg_signal(FATAL,"The 'post' action expects: chidg post gridfile.h5 solutionfile.h5")
                call get_command_argument(2,grid_file)
                call get_command_argument(3,solution_file)
                call chidg_post(trim(grid_file), trim(solution_file))
                call chidg_post_vtk(trim(grid_file), trim(solution_file))

            case ('matplotlib')
                if (narg /= 3) call chidg_signal(FATAL,"The 'matplotlib' action expects: chidg matplotlib gridfile.h5solutionfile.h5")
                call get_command_argument(2,grid_file)
                call get_command_argument(3,solution_file)
                call chidg_post_matplotlib(trim(grid_file),trim(solution_file))

            case ('airfoil')
                if (narg /= 2) call chidg_signal(FATAL,"The 'airfoil' action expects: chidg airfoil solutionfile.h5")
                call get_command_argument(2,solution_file)
                call chidg_airfoil(trim(solution_file))

            case ('clone')
                if (narg /= 3) call chidg_signal(FATAL,"The 'clone' action expects: chidg clone source_file.h5 target_file.h5")
                call get_command_argument(2,file_a)
                call get_command_argument(3,file_b)
                call chidg_clone(trim(file_a),trim(file_b))

            case default
                call chidg_signal(FATAL,"We didn't understand the way chidg was called. Available chidg 'actions' are: 'edit' 'convert' 'post' 'matplotlib' and 'airfoil'.")
        end select

        call chidg%shut_down('core')





    else
        call chidg_signal(FATAL,"chidg: invalid number of arguments. Expecting (0) arguments: 'chidg'. or (2) arguments: 'chidg action file'.")
    end if












end program driver
