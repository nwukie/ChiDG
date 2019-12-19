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
    use mod_constants,          only: OUTPUT_RES
    use type_chidg,             only: chidg_t
    use type_dict,              only: dict_t
    use mod_tecio,              only: write_tecio_file
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
    subroutine chidg_post_hdf2tec(chidg,grid_file,solution_file)
        type(chidg_t),  intent(inout)   :: chidg
        character(*),   intent(in)      :: grid_file
        character(*),   intent(in)      :: solution_file
    
        type(file_properties_t)             :: file_props
        character(:),           allocatable :: time_string, solution_file_prefix, plt_filename
        integer(ik)                         :: nterms_s, solution_order, ierr, iproc

        ! Get nterms_s 
        do iproc = 0,NRANK-1
            if (iproc == IRANK) then
                file_props  = get_properties_hdf(solution_file)
                nterms_s    = file_props%nterms_s(1)
                time_string = file_props%time_integrator

                solution_order = 0
                do while (solution_order*solution_order*solution_order < nterms_s)
                    solution_order = solution_order + 1
                end do
            end if
            call MPI_Barrier(ChiDG_COMM,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'chidg_post_hdf2tec_new: error in MPI_Barrier.')
        end do


        ! Initialize solution data storage
        call chidg%set('Solution Order', integer_input=solution_order)
        call chidg%set('Time Integrator', algorithm=trim(time_string))
        chidg%grid_file        = grid_file
        chidg%solution_file_in = solution_file


        ! Read grid/solution modes and time integrator options from HDF5
        !call chidg%read_mesh(grid_file, 'primal storage', interpolation='Uniform', level=OUTPUT_RES)
        call chidg%read_mesh(grid_file, 'adjoint storage', interpolation='Uniform', level=OUTPUT_RES)


        call chidg%data%set_fields_post(file_props%contains_primary_fields,file_props%contains_adjoint_fields,file_props%nfunctionals)

        if (file_props%contains_primary_fields) call chidg%read_fields(solution_file,'primary')
        if (file_props%contains_adjoint_fields) call chidg%read_fields(solution_file,'adjoint')


        ! Process for getting wall distance
        call chidg%process()


        ! Get post processing data (q_out)
        call chidg%time_integrator%initialize_state(chidg%data)
        do iproc = 0,NRANK-1
            if (iproc == IRANK) then
                call chidg%time_integrator%read_time_options(chidg%data,solution_file,'process')
            end if
            call MPI_Barrier(ChiDG_COMM,ierr)
            if (ierr /= 0) call chidg_signal(FATAL,'chidg_post_hdf2tec_new: error in MPI_Barrier.')
        end do

        ! Move solution vector from q_in to q
        if (file_props%contains_primary_fields) call chidg%time_integrator%process_data_for_output(chidg%data)
        ! Move solution vector from v_in to v, v is initialized here for post-processing
        if (file_props%contains_adjoint_fields) call chidg%data%sdata%adjoint%process_adjoint_solution()


        ! Write solution
        solution_file_prefix = get_file_prefix(solution_file,'.h5')
        call write_tecio_file(chidg%data,solution_file_prefix, write_domains=.true., write_surfaces=.true.)
        

    end subroutine chidg_post_hdf2tec
    !******************************************************************************************
    



end module mod_chidg_post_hdf2tec
