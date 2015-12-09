!> Chimera-based, discontinuous Galerkin equation solver
!!
!! This program is designed to solve partial differential equations,
!! and systems of partial differential equations, using the discontinuous
!! Galerkin method for spatial discretization using Chimera, overset grids to
!! represent the simulation domain.
!!
!! @author Nathan A. Wukie
!!
!!
!---------------------------------------------------------------------------------------------



program driver
    use mod_kinds,              only: rk, ik
    use type_chidg,             only: chidg_t
    use type_dict,              only: dict_t
    use mod_tecio,              only: write_tecio_variables
    
    !
    ! Variable declarations
    !
    implicit none
    type(chidg_t)                       :: chidg

    character(len=:), allocatable       :: gridfile, solutionfile
    integer                             :: nargs



    !
    ! Initialize ChiDG environment
    !
    print*, 'chidg init'
    call chidg%init('env')
    !call chidg%init('io')



    !
    ! Get number of arguments
    !
    nargs = command_argument_count()



    print*, 'hi  - 2'
    !
    ! Check if a filename was provided to the program
    !
    if ( nargs < 2) then
        print*, "Usage: H5toTEC 'gridfile' 'solutionfile'"
        stop
    end if



    print*, 'hi  - 3'
    !
    ! Get file name from command-line argument
    !
    call get_command_argument(1, gridfile)
    call get_command_argument(2, solutionfile)





    print*, 'hi  - 4'


    !
    ! Get nterms_s and eqnset
    !
    solution_properties = get_file_properties(solutionfile)

    nterms_s    = solution_properties%nterms_s
    equationset = solution_properties%equationset


    !
    ! Read grid data from file
    !
    call chidg%read_grid(gridfile)






    print*, 'hi  - 5'
    !
    ! Initialize solution data storage
    !
    call chidg%init('chimera')
    call chidg%data%init_sdata()




    print*, 'hi  - 6'
    !
    ! Read solution modes from HDF5
    !
    call chidg%read_solution(solutionfile)




    print*, 'hi  - 7'
    !
    ! Write solution in TecIO format
    !
    call write_tecio_variables(chidg%data,'0.plt',1)


    

    !
    ! Close ChiDG
    !
    call chidg%close()





end program driver
