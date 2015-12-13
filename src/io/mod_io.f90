module mod_io
#include <messenger.h>
    !!          Module Includes             !!
    use mod_kinds,     only: rk,ik
    use mod_constants, only: MAXBLOCKS

    !!          Variable Declarations       !!
    implicit none


    ! FILES
    !--------------------------------------------------
    character(len=100),  save    :: gridfile
    character(len=100),  save    :: gridtype
    character(len=100),  save    :: tecplot_prefix      = 'tec'
    character(len=100),  save    :: hdf_out             = 'solution.h5'

    character(len=100),  save    :: solutionfile_in     = 'none'
    character(len=100),  save    :: solutionfile_out    = 'none'




    ! SPACE
    !--------------------------------------------------
    character(len=100),  save    :: basis = 'legendre'
    integer(ik),         save    :: solution_order  = 1

 


    
    ! QUADRATURE
    !--------------------------------------------------
    integer(ik),         save    :: gq_rule = 2          !> 1: Collocation, 2: Over-integration
   


   
    
    ! EQUATION SET 
    !--------------------------------------------------
    character(len=100),  save    :: eqnset = 'scalar'
  
  



    
    ! TIME
    !--------------------------------------------------
    character(len=100),  save    :: timescheme
    real(rk),            save    :: dt = 0.001_rk
    real(rk),            save    :: cfl0 = 1._rk
    integer(ik),         save    :: nsteps = 100
    real(rk),            save    :: ttol = 1.e-8
    integer(ik),         save    :: ntime_instances = 1
   
   


   
    ! MATRIX SOLVER
    !--------------------------------------------------
    character(len=100),  save    :: matrixsolver = 'direct'
    real(rk),            save    :: mtol = 1.e-8



    ! PRECONDITIONER
    !--------------------------------------------------
    character(len=100),  save    :: preconditioner = 'identity'
   
   



    ! IO
    !--------------------------------------------------
    integer(ik),         save    :: nwrite         = 100
    logical,             save    :: initial_write  = .true.
    logical,             save    :: final_write    = .true.
    integer(ik),         save    :: output_res     = 5
     
    



    !==================================================================================
    !           These quantities are used globally, but computed during input.
    !           So, they do not need explicitly initialized in the namelist file
    !==================================================================================
    integer(ik),         save    :: nterms_sol1d    = 1
    integer(ik),         save    :: nterms_s        = 1

contains






    !> Read input file
    !!      - Read chidg.nml namelist formatted input file
    !!
    !!
    !!  @author Nathan A. Wukie
    !!
    !------------------------------------------------------------------------------------------------
    subroutine read_input()
        use mod_ordering,   only: MAX_POLY_ORDER

        logical :: file_exists

        namelist /files/                    gridfile,              &
                                            gridtype,              &
                                            hdf_out,               &
                                            tecplot_prefix,        &
                                            solutionfile_in,       &
                                            solutionfile_out

        namelist /space/                    basis,                 &
                                            solution_order

        namelist /quadrature/               gq_rule

        namelist /equation_set/             eqnset

        namelist /time/                     timescheme,            &
                                            cfl0,                  &
                                            dt,                    &
                                            nsteps,                &
                                            ntime_instances,       &
                                            ttol

        namelist /matrix_solver/            matrixsolver,          &
                                            mtol,                  &
                                            preconditioner


        namelist /io/                       nwrite,                &
                                            output_res,            &
                                            initial_write,         &
                                            final_write


        !
        ! Check that input file exists
        !
        inquire(file='chidg.nml', exist=file_exists)
        if (.not. file_exists) call chidg_signal(FATAL, "read_input: 'chidg.nml' input file was not found")



        !
        ! Read namelist input for parameter initialization
        !
        open(unit=7,form='formatted',file="chidg.nml")
        read(7,nml=files)
        read(7,nml=space)
        read(7,nml=quadrature)
        read(7,nml=equation_set)
        read(7,nml=time)
        read(7,nml=matrix_solver)
        read(7,nml=io)



        !
        ! Compute number of terms in polynomial expansions
        !
        nterms_sol1d = (solution_order)
        nterms_s = nterms_sol1d * nterms_sol1d * nterms_sol1d


    end subroutine read_input
    !#######################################################################################################




end module mod_io
