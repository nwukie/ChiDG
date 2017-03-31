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
        integer(ik)                         :: nterms_s, spacedim, solution_order



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




        !
        ! Read grid data from file
        !
        call chidg%read_grid(grid_file,spacedim)


        solution_order = 0
        do while (solution_order*solution_order*solution_order < nterms_s)
            solution_order = solution_order + 1
        end do


        !
        ! Initialize solution data storage
        !
        call chidg%set('Solution Order', integer_input=solution_order)
        call chidg%init('domains')
        call chidg%init('communication')
        call chidg%init('solvers')




        !
        ! Read solution modes from HDF5
        !
        call chidg%read_solution(solution_file)




        !
        ! Write solution in TecIO format
        !
        call write_tecio_variables(chidg%data,'0.plt')


        

        !
        ! Close ChiDG
        !
        call chidg%shut_down('core')



    end subroutine chidg_post_hdf2tec
    !******************************************************************************************


    



end module mod_chidg_post_hdf2tec
