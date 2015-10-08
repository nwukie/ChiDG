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
    character(len=100),  save    :: tecplot_prefix = 'tec'
    character(len=100),  save    :: hdf_out        = 'solution.h5'

    character(len=100),  save    :: solutionfile




    ! SPACE
    !--------------------------------------------------
    character(len=100),  save    :: basis = 'legendre'
    integer(ik),         save    :: solution_order  = 1

 


    
    ! QUADRATURE
    !--------------------------------------------------
    integer(ik),         save    :: gq_rule = 2          !> 1: Collocation, 2: Over-integration
   


   
    
    ! EQUATION SET 
    !--------------------------------------------------
    character(len=100),  save    :: eqnset
  
  



    
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
    integer(ik),         save    :: output_res     = 5
     
    


    ! BOUNDARY CONDITIONS
    !--------------------------------------------------
    integer(ik),         save    :: bc_ximin(MAXBLOCKS),   bc_ximax(MAXBLOCKS)
    integer(ik),         save    :: bc_etamin(MAXBLOCKS),  bc_etamax(MAXBLOCKS)
    integer(ik),         save    :: bc_zetamin(MAXBLOCKS), bc_zetamax(MAXBLOCKS)
    real(rk),            save    :: bcpar1(6), bcpar2(6), bcpar3(6), bcpar4(6)


    !==================================================================================
    !           These quantities are used globally, but computed during input.
    !           So, they do not need explicitly initialized in the namelist file
    !==================================================================================
    integer(ik),         save    :: nterms_sol1d    = 1
    integer(ik),         save    :: nterms_s        = 1

contains
!--------------------------------------------------






    !> Read input file
    !!      - Read chidg.nml namelist formatted input file
    !!
    !!
    !!  @author Nathan A. Wukie
    !!
    !--------------------------------------------------------
    subroutine read_input()
        use mod_ordering,   only: MAX_POLY_ORDER

        logical :: file_exists

        namelist /files/                    gridfile,              &
                                            gridtype,              &
                                            hdf_out,               &
                                            tecplot_prefix


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

!        namelist /grid/     gridfile,   gridtype,   &
!                            bc_ximin,   bc_ximax,   &
!                            bc_etamin,  bc_etamax,  &
!                            bc_zetamin, bc_zetamax, &
!                            bcpar1, bcpar2, bcpar3, bcpar4

        namelist /io/       nwrite, output_res

        inquire(file='chidg.nml', exist=file_exists)
        if (.not. file_exists) call signal(FATAL, "read_input: 'chidg.nml' input file was not found")

        open(unit=7,form='formatted',file="chidg.nml")
        read(7,nml=files)
        read(7,nml=space)
        read(7,nml=quadrature)
        read(7,nml=equation_set)
        read(7,nml=time)
        read(7,nml=matrix_solver)
        read(7,nml=io)



        ! Compute number of terms in polynomial expansions
        nterms_sol1d = (solution_order)
        nterms_s = nterms_sol1d * nterms_sol1d * nterms_sol1d


    end subroutine read_input




end module mod_io
