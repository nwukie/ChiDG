!> ChiDG HDF5 to Tecplot TecIO conversion utility.
!!
!!  Given HDF5 files containing ChiDG formatted grid and solution data, H5toTEC will produce
!!  a tecplot file for visualizing the solution by sampling the solution and coordinate polynomials.
!!
!!  The tecplot file produced cannot be used in any way other than for visualization purposes.
!!
!! @author Nathan A. Wukie
!!
!!
!! Usage:   H5toTEC 'gridfile' 'solutionfile'
!!
!!
!---------------------------------------------------------------------------------------------
program HDF5toTEC
#include <messenger.h>
    use mod_kinds,              only: rk, ik
    use type_chidg,             only: chidg_t
    use type_dict,              only: dict_t
    use mod_tecio,              only: write_tecio_variables
    use type_file_properties,   only: file_properties_t
    use mod_file_utilities,     only: get_file_properties
    use mod_io,                 only: nterms_s, eqnset
    
    !
    ! Variable declarations
    !
    implicit none
    type(chidg_t)                       :: chidg

    character(len=1024)                 :: gridfile, solutionfile
    integer                             :: nargs

    type(file_properties_t)             :: file_props



    !
    ! Initialize ChiDG environment
    !
    call chidg%init('env')
    !call chidg%init('io')




    !
    ! Get number of arguments
    !
    nargs = command_argument_count()




    !
    ! Check if a filename was provided to the program
    !
    if ( nargs < 2) call chidg_signal(FATAL,"Usage: H5toTEC 'gridfile' 'solutionfile'")




    !
    ! Get file name from command-line argument
    !
    call get_command_argument(1, gridfile)
    call get_command_argument(2, solutionfile)





    !
    ! Get nterms_s and eqnset. TODO: Read directly from file
    !
    file_props = get_file_properties(solutionfile)

    nterms_s    = file_props%nterms_s(1)    ! Global variable from mod_io
    eqnset      = file_props%eqnset(1)      ! Global variable from mod_io




    !
    ! Read grid data from file
    !
    call chidg%read_grid(gridfile)




    !
    ! Initialize solution data storage
    !
!    call chidg%init('chimera')
    call chidg%data%init_sdata()




    !
    ! Read solution modes from HDF5
    !
    call chidg%read_solution(solutionfile)




    !
    ! Write solution in TecIO format
    !
    call write_tecio_variables(chidg%data,'0.plt',1)


    

    !
    ! Close ChiDG
    !
    call chidg%close()


    



end program HDF5toTEC
