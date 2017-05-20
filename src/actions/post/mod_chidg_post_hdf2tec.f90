!> ChiDG HDF5 to Tecplot TecIO conversion utility.
!!
!!  Given HDF5 files containing ChiDG formatted grid and solution data, HDF5toTEC will produce
!!  a tecplot file for visualizing the solution by sampling the solution and coordinate polynomials.
!!
!!  The tecplot file produced cannot be used in any way other than for visualization purposes.
!!
!!  @author Nathan A. Wukie
!!  @date   2/3/2016
!!
!!
!! Usage:   chidg post 'chidgfile'
!!
!!
!---------------------------------------------------------------------------------------------
module mod_chidg_post_hdf2tec
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use type_chidg,             only: chidg_t
    use type_dict,              only: dict_t
    use mod_tecio,              only: write_tecio_variables
    use type_file_properties,   only: file_properties_t
    use mod_hdf_utilities,      only: get_properties_hdf
    use mod_string,             only: get_file_prefix
    implicit none







contains



    !>  Post-processing tool for writing a tecplot file of sampled modal data.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/6/2016
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------
    subroutine chidg_post_hdf2tec(grid_file,solution_file)
        character(*)    :: grid_file
        character(*)    :: solution_file
    
        type(chidg_t)                       :: chidg
        type(file_properties_t)             :: file_props
        character(:),           allocatable :: eqnset
        character(:),           allocatable :: time_string, grid_file_prefix, solution_file_prefix, step_str, &
                                               plt_filename
        integer(ik)                         :: nterms_s, spacedim, solution_order, istep, itimestep



        !
        ! Initialize ChiDG environment
        !
        call chidg%start_up('core')
        call chidg%start_up('mpi')




        !
        ! Get nterms_s and eqnset.
        !
        file_props = get_properties_hdf(solution_file)

        nterms_s    = file_props%nterms_s(1)
        eqnset      = file_props%eqnset(1)
        spacedim    = file_props%spacedim(1)
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
        step_str             = solution_file_prefix(len(grid_file_prefix)+2:)

        plt_filename = step_str//'.plt'

        !
        ! Write solution
        !
        call write_tecio_variables(chidg%data,plt_filename)
        

        !
        ! Close ChiDG
        !
        call chidg%shut_down('core')



    end subroutine chidg_post_hdf2tec
    !******************************************************************************************


    



end module mod_chidg_post_hdf2tec
