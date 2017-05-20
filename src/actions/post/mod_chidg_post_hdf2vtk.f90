!>  ChiDG HDF5 to vtk conversion utility
!!
!!  @author Mayank Sharma
!!  @date 10/31/2016
!!
!!
!!  Usage: chidg post 'chidgfile'
!!  
!
!----------------------------------------------------------------------------------------------
module mod_chidg_post_hdf2vtk
#include <messenger.h>
    use mod_kinds,              only: rk,ik
    use type_chidg,             only: chidg_t
    use type_dict,              only: dict_t
    use mod_vtkio,              only: write_vtk_file
    use type_file_properties,   only: file_properties_t
    use mod_hdf_utilities,      only: get_properties_hdf
    use mod_string,             only: get_file_prefix

    implicit none



contains



    !>  Post processing tool for writing a vtk file of sampled modal data
    !!
    !!  @author Mayank Sharma
    !!  @date 10/31/2016
    !!
    !!  Functionality for unsteady IO
    !!
    !!  @author Mayank Sharma
    !!  @date   3/22/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_post_hdf2vtk(grid_file,solution_file)
        character(*)                    ::  grid_file
        character(*)                    ::  solution_file


        type(chidg_t)                       ::  chidg
        type(file_properties_t)             ::  file_props
        character(:),           allocatable ::  eqnset
        character(:),           allocatable ::  time_string, grid_file_prefix, solution_file_prefix, step_str, &
                                                pvd_filename, time_integrator_string
        integer(ik)                         ::  nterms_s,spacedim,solution_order,istep,itimestep


        !
        ! Initialize ChiDG environment
        !
        call chidg%start_up('mpi')
        call chidg%start_up('core')



        !
        ! Get nterms_s and eqnset
        !
        file_props = get_properties_hdf(solution_file)

        nterms_s    = file_props%nterms_s(1)     ! Global variable from mod_io 
        eqnset      = file_props%eqnset(1)       ! Global variable from mod_io
        spacedim    = file_props%spacedim(1)     ! Global variable from mod_io
        time_string = file_props%time_integrator






        solution_order = 0
        do while (solution_order*solution_order*solution_order < nterms_s)
            solution_order = solution_order + 1
        end do


        !
        ! Initialize solution data storage
        !
        call chidg%set('Solution Order', integer_input=solution_order)
        call chidg%set('Time Integrator', algorithm=trim(time_string))
        call chidg%time_integrator%initialize_state(chidg%data)
        !call chidg%init('domains')
        !call chidg%init('communication')
        !call chidg%init('solvers')


        !
        ! Read grid/solution modes and time integrator options from HDF5
        !
        call chidg%read_grid(grid_file,spacedim)
        call chidg%read_solution(solution_file)
        call chidg%time_integrator%read_time_options(chidg%data,solution_file)

        
        !
        ! Get post processing data (q_out)
        !        
        call chidg%time_integrator%process_data_for_output(chidg%data)
       
        
        grid_file_prefix     = get_file_prefix(grid_file,'.h5')
        solution_file_prefix = get_file_prefix(solution_file,'.h5')
        step_str = solution_file_prefix(len(grid_file_prefix) + 2:)
        read(step_str,*) istep
        itimestep = istep/chidg%data%time_manager%nwrite
       
        pvd_filename = 'chidg_results.pvd'
        

        !
        ! Write solution in vtk format
        !
        call write_vtk_file(chidg%data,itimestep,pvd_filename)


        !
        ! Shut down ChiDG
        !
        call chidg%shut_down('core')


    end subroutine chidg_post_hdf2vtk 
    !******************************************************************************************




















end module mod_chidg_post_hdf2vtk
