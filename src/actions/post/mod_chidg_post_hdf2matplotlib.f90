!>  ChiDG HDF5 to matplotlib conversion utility
!!
!!  @author Mayank Sharma
!!  @date   3/24/2017
!!
!!  Usage: chidg post 'chidgfile'
!!
!-------------------------------------------------------------------------------------------------
module mod_chidg_post_hdf2matplotlib
#include<messenger.h>
    use mod_kinds,              only: rk,ik
    use type_chidg,             only: chidg_t
    use type_dict,              only: dict_t
    use mod_matplotlib_io,      only: write_matplotlib_file
    use type_file_properties,   only: file_properties_t
    use mod_hdf_utilities,      only: get_properties_hdf

    implicit none



contains



    !>  Post processing tool for writing a .dat file of sampled modal data for matplotlib
    !!
    !!  @author Mayank Sharma
    !!  @date   3/24/2017
    !!
    !---------------------------------------------------------------------------------------------
    subroutine chidg_post_hdf2matplotlib(filename)
        character(*)                :: filename

        type(chidg_t)               :: chidg
        type(file_properties_t)     :: file_props
        character(:),   allocatable :: eqnset
        character(:),   allocatable :: time_string, coord, matplotlib_file, &
                                       original_sol_file, Fourier_coeff_file
        integer(ik)                 :: nterms_s, spacedim, solution_order


        !
        ! Initialize ChiDG environment
        !
        call chidg%start_up('mpi')
        call chidg%start_up('core')


        !
        ! Get nterms_s and eqnset
        !
        file_props = get_properties_hdf(filename)

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
        call chidg%read_grid(filename,spacedim)
        call chidg%read_solution(filename)
        call chidg%time_integrator%read_time_options(chidg%data,filename)


        !
        ! Get post processing data (q_out)
        !
        call chidg%time_integrator%process_data_for_output(chidg%data)


        !
        ! Write solution for matplotlib
        !
        coord              = 'x'
        matplotlib_file    = 'chidg_results.dat'
        original_sol_file  = 'chidg_original_results.dat'
        Fourier_coeff_file = 'chidg_coefficients.dat'
        call write_matplotlib_file(chidg%data,1,coord,matplotlib_file,filename_2 = original_sol_file,&
                                   filename_3 = Fourier_coeff_file)


        !
        ! Shut down ChiDG
        !
        call chidg%shut_down('core')


    end subroutine chidg_post_hdf2matplotlib
    !*********************************************************************************************




















end module mod_chidg_post_hdf2matplotlib
