!>  Chimera-based, discontinuous Galerkin equation solver
!!
!!  This program is designed to solve partial differential equations,
!!  and systems of partial differential equations, using the discontinuous
!!  Galerkin method for spatial discretization using overset grids to
!!  represent the simulation domain.
!!
!!  @author Nathan A. Wukie
!!  @date   1/31/2016
!!
!---------------------------------------------------------------------------------------------
program driver
#include <messenger.h>
    use type_chidg,                 only: chidg_t
    use type_function,              only: function_t
    use mod_function,               only: create_function
    use mod_file_utilities,         only: delete_file
    use mod_io

    ! Actions
    use mod_chidg_edit,         only: chidg_edit
    use mod_chidg_convert,      only: chidg_convert
    use mod_chidg_post,         only: chidg_post, chidg_post_vtk, chidg_post_matplotlib
    use mod_chidg_forces,       only: chidg_forces
    use mod_chidg_clone,        only: chidg_clone
    use mod_chidg_post_hdf2tec, only: chidg_post_hdf2tec_new
    use mod_tutorials,          only: tutorial_driver
    use mod_euler_eigenmodes,   only: compute_euler_eigenmodes

    ! Variable declarations
    implicit none
    type(chidg_t)                               :: chidg
    integer                                     :: narg, ierr, ifield
    character(len=1024)                         :: chidg_action, filename, grid_file, solution_file, file_a, file_b, file_in, pattern, tutorial, patch_group
    character(len=10)                           :: time_string
    character(:),                   allocatable :: command, tmp_file
    class(function_t),              allocatable :: fcn


    ! Check for command-line arguments
    narg = command_argument_count()


    ! Execute ChiDG calculation
    if ( narg == 0 ) then

        ! Initialize ChiDG environment
        call chidg%start_up('mpi')
        call chidg%start_up('namelist')
        call chidg%start_up('core')


        ! Set ChiDG Algorithms, Accuracy
        call chidg%set('Time Integrator' , algorithm=time_integrator                   )
        call chidg%set('Nonlinear Solver', algorithm=nonlinear_solver, options=noptions)
        call chidg%set('Linear Solver'   , algorithm=linear_solver,    options=loptions)
        call chidg%set('Preconditioner'  , algorithm=preconditioner                    )
        call chidg%set('Solution Order'  , integer_input=solution_order                )


        ! Read grid and boundary condition data
        call chidg%read_mesh(gridfile)


        ! Initialize solution
        !   1: 'none', init fields with values from mod_io module variable initial_fields(:)
        !   2: read initial solution from ChiDG hdf5 file
        if (solutionfile_in == 'none') then
            call create_function(fcn,'constant')
            do ifield = 1,chidg%data%mesh%domain(1)%neqns
                call fcn%set_option('val',initial_fields(ifield))
                call chidg%data%sdata%q_in%project(chidg%data%mesh,fcn,ifield)
            end do

        else
            call chidg%read_fields(solutionfile_in)
        end if


        ! Run ChiDG simulation
        call chidg%report('before')
        call chidg%run(write_initial=initial_write, write_final=final_write)
        call chidg%report('after')


        ! Close ChiDG
        call chidg%shut_down('core')
        call chidg%shut_down('mpi')





    ! Check if executing 'action'
    else if ( narg >= 1 ) then

        ! Get 'action'
        call get_command_argument(1,chidg_action)
        !call chidg%start_up('core')

        ! Select 'action'
        select case (trim(chidg_action))
            !>  ChiDG:convert   src/actions/convert
            !!
            !!  Convert Multi-block, Unformatted, Double-Precision, Plot3D grids to
            !!  ChiDG-formatted HDF5 file.
            !!
            !!  NOTE: this routine handles agglomeration of linear elements to form 
            !!  higher-order elements.
            !!
            !!  Command-Line:
            !!  --------------------
            !!  chidg convert myfile.x
            !!
            !!  Produces:
            !!  --------------------
            !!  myfile.h5
            !!
            !----------------------------------------------------------------------------
            case ('convert')
                call chidg%start_up('core')
                if (narg /= 2) call chidg_signal(FATAL,"The 'convert' action expects: chidg convert filename.x")
                call get_command_argument(2,filename)
                call chidg_convert(trim(filename))

            !*****************************************************************************



            !>  ChiDG:edit  src/actions/edit
            !!
            !!  Edit a ChiDG HDF5 file. Edit equations, boundary conditions + settings,
            !!  and patches.
            !!
            !!  Command-Line:
            !!  ---------------------
            !!  chidg edit myfile.h5
            !!
            !----------------------------------------------------------------------------
            case ('edit')
                call chidg%start_up('core')
                call get_command_argument(2,filename)
                if (narg < 2) call chidg_signal(FATAL,"The 'edit' action was called with too few arguments. Try: chidg edit filename.h5")
                if (narg > 2) call chidg_signal(FATAL,"The 'edit' action was called with too many arguments. Try: chidg edit filename.h5")
                call chidg_edit(trim(filename))

            !*****************************************************************************



            case ('post')
            !>  ChiDG:post  src/actions/post
            !!
            !!  Post-process solution files for visualization.
            !!
            !!  Command-Line MODE 1: Single-file
            !!  --------------------------------
            !!
            !!     Command-line:                    Output:
            !!  chidg post myfile.h5       myfile.plt (Tecplot-readable)
            !!
            !!  Command-Line MODE 2: Multi-file
            !!  --------------------------------
            !!  In the case where there are several files that need processed,
            !!  wildcards can be passed in, but must be wrapped in quotes " ".
            !!
            !!  Files: myfile_0.1000.h5, myfile_0.2000.h5, myfile_0.3000.h5
            !!
            !!     Command-line:                Output:
            !!  chidg post "myfile*"        myfile_0.1000.plt
            !!                              myfile_0.2000.plt
            !!                              myfile_0.3000.plt
            !!
            !!---------------------------------------------------------------------------
                call chidg%start_up('core')
                if (narg /= 2) call chidg_signal(FATAL,"The 'post' action expects: chidg post file.h5")

                call date_and_time(time=time_string)
                tmp_file = 'chidg_post_files'//time_string//'.txt'
                call get_command_argument(2,pattern)
                command = 'ls '//trim(pattern)//' > '//tmp_file
                call system(command)
            

                open(7,file=tmp_file,action='read')
                do
                    read(7,fmt='(a)', iostat=ierr) solution_file
                    if (ierr /= 0) exit
                    call chidg_post(trim(solution_file), trim(solution_file))
                end do
                close(7)

                call delete_file(tmp_file)
            !*****************************************************************************



            case ('2tec')
            !>  ChiDG:post  src/actions/post
            !!
            !!  Post-process solution files for visualization (tecplot format)
            !!
            !!  Command-Line MODE 1: Single-file
            !!  --------------------------------
            !!
            !!     Command-line:                    Output:
            !!  chidg post myfile.h5       myfile.plt (Tecplot-readable)
            !!
            !!  Command-Line MODE 2: Multi-file
            !!  --------------------------------
            !!  In the case where there are several files that need processed,
            !!  wildcards can be passed in, but must be wrapped in quotes " ".
            !!
            !!  Files: myfile_0.1000.h5, myfile_0.2000.h5, myfile_0.3000.h5
            !!
            !!     Command-line:                Output:
            !!  chidg post "myfile*"        myfile_0.1000.plt
            !!                              myfile_0.2000.plt
            !!                              myfile_0.3000.plt
            !!
            !!---------------------------------------------------------------------------
                call chidg%start_up('core')
                if (narg /= 2) call chidg_signal(FATAL,"The '2tec' action expects: chidg 2tec file.h5")

                call date_and_time(time=time_string)
                tmp_file = 'chidg_2tec_files'//time_string//'.txt'
                call get_command_argument(2,pattern)
                command = 'ls '//trim(pattern)//' > '//tmp_file
                call system(command)
            

                open(7,file=tmp_file,action='read')
                do
                    read(7,fmt='(a)', iostat=ierr) solution_file
                    if (ierr /= 0) exit
                    call chidg_post_hdf2tec_new(trim(solution_file),trim(solution_file))
                end do
                close(7)

                call delete_file(tmp_file)
            !*****************************************************************************



            case ('2vtk')
            !>  ChiDG:post  src/actions/post
            !!
            !!  Post-process solution files for visualization (vtk format).
            !!
            !!  Command-Line MODE 1: Single-file
            !!  --------------------------------
            !!
            !!     Command-line:                    Output:
            !!  chidg post myfile.h5       myfile_itime_idom_itimestep.vtu (Paraview-readable)
            !!
            !!  Command-Line MODE 2: Multi-file
            !!  --------------------------------
            !!  In the case where there are several files that need processed,
            !!  wildcards can be passed in, but must be wrapped in quotes " ".
            !!
            !!  Files: myfile_0.1000.h5, myfile_0.2000.h5, myfile_0.3000.h5
            !!
            !!     Command-line:                Output:
            !!  chidg post "myfile*"        myfile_0_0_1.vtu
            !!                              myfile_0_0_2.vtu
            !!                              myfile_0_0_3.vtu
            !!
            !!---------------------------------------------------------------------------
                call chidg%start_up('core')
                if (narg /= 2) call chidg_signal(FATAL,"The '2vtk' action expects: chidg 2vtk file.h5")

                call date_and_time(time=time_string)
                tmp_file = 'chidg_2vtk_files'//time_string//'.txt'
                call get_command_argument(2,pattern)
                command = 'ls '//trim(pattern)//' > '//tmp_file
                call system(command)
            

                open(7,file=tmp_file,action='read')
                do
                    read(7,fmt='(a)', iostat=ierr) solution_file
                    if (ierr /= 0) exit
                    call chidg_post_vtk(trim(solution_file), trim(solution_file))
                end do
                close(7)

                call delete_file(tmp_file)
            !*****************************************************************************


    
            !>  ChiDG:clone src/actions/clone
            !!
            !!  Clone a ChiDG-file configuration from one file to another.
            !!
            !!  Command-Line:
            !!  ------------------------
            !!  chidg clone source.h5 target.h5
            !!
            !!  MODE1: Copy boundary condition state groups AND patch attributes 
            !!         (assumes the grid domain/topology/names match from source to target.
            !!  MODE2: Copy boundary condition state groups ONLY
            !!  MODE3: Copy patch attributes ONLY
            !!         (assumes the grid domain/topology/names match from source to target.
            !!
            !-----------------------------------------------------------------------------
            case ('clone')
                call chidg%start_up('core')
                if (narg /= 3) call chidg_signal(FATAL,"The 'clone' action expects: chidg clone source_file.h5 target_file.h5")
                call get_command_argument(2,file_a)
                call get_command_argument(3,file_b)
                call chidg_clone(trim(file_a),trim(file_b))

            !*****************************************************************************




            case ('matplotlib')
                call chidg%start_up('core')
                if (narg /= 3) call chidg_signal(FATAL,"The 'matplotlib' action expects: chidg matplotlib gridfile.h5solutionfile.h5")
                call get_command_argument(2,grid_file)
                call get_command_argument(3,solution_file)
                call chidg_post_matplotlib(trim(grid_file),trim(solution_file))

            case ('forces')
                call chidg%start_up('core')
                if (narg /= 2) call chidg_signal(FATAL,"The 'forces' action expects to be called as: chidg forces solutionfile.h5")
                call get_command_argument(2,solution_file)
                call write_line('Enter patch group to integrate: ')
                read*, patch_group


                call date_and_time(time=time_string)
                tmp_file = 'chidg_forces_files'//time_string//'.txt'
                call get_command_argument(2,pattern)
                command = 'ls '//trim(pattern)//' > '//tmp_file
                call system(command)
            
                open(7,file=tmp_file,action='read')
                do
                    read(7,fmt='(a)', iostat=ierr) solution_file
                    if (ierr /= 0) exit
                    call chidg_forces(trim(solution_file),trim(patch_group))
                end do
                close(7)

                call delete_file(tmp_file)



            case ('inputs')
                call chidg%start_up('core')
                if (narg > 2) call chidg_signal(FATAL,"The 'inputs' action expects to be called as: chidg inputs")
                call write_namelist()


            case ('tutorial')
                call chidg%start_up('core')
                if (narg /= 2) call chidg_signal(FATAL,"The 'tutorial' action expects to be called as: chidg tutorial selected_tutorial.")
                call get_command_argument(2,tutorial)
                call tutorial_driver(trim(tutorial))

            case ('eigen')
                call chidg%start_up('core')
                call compute_euler_eigenmodes()

            case default
                call chidg%start_up('core')
                call chidg_signal(FATAL,"We didn't understand the way chidg was called. Available chidg 'actions' are: 'edit' 'convert' 'post' 'matplotlib' 'inputs' and 'forces'.")
        end select

        call chidg%shut_down('core')





    else
        call chidg_signal(FATAL,"chidg: invalid number of arguments. Expecting (0) arguments: 'chidg'. or (2) arguments: 'chidg action file'.")
    end if











end program driver
