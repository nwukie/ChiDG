module mod_io
#include <messenger.h>
    !!          Module Includes             !!
    use mod_kinds,     only: rk,ik
    use mod_constants, only: MAXBLOCKS

    !!          Variable Declarations       !!
    implicit none


    ! POLYNOMIAL
    !--------------------------------------------------
    character(len=100),  save    :: basis = 'legendre'


    ! GRID  
    !--------------------------------------------------
    character(len=100),  save    :: gridfile
    character(len=100),  save    :: gridtype


    ! SOLUTION
    !--------------------------------------------------
    character(len=100),  save    :: solutionfile
    integer(ik),         save    :: solution_order  = 1
    integer(ik),         save    :: nterms_sol1d    = 1
    integer(ik),         save    :: nterms_s        = 1
    integer(ik),         save    :: ntime_instances = 1

    !integer(ik),         save    :: nterms_sol1d = 2
    !integer(ik),         save    :: nterms_sol2d = 2
    !integer(ik),         save    :: nterms_sol3d = 2
    !integer(ik),         save    :: nterms_mesh1d = 1
    !integer(ik),         save    :: nterms_mesh2d = 4
    !integer(ik),         save    :: nterms_mesh3d = 8
 
    
    ! QUADRATURE
    !--------------------------------------------------
    integer(ik),         save    :: gq_rule = 2          !> 1: Collocation, 2: Over-integration
   
   
    
    ! EQUATION SET 
    !--------------------------------------------------
    character(len=100),  save    :: eqnset
  
  
    
    ! SOLUTION ADVANCEMENT 
    !--------------------------------------------------
    character(len=100),  save    :: temporal_scheme
    real(rk),            save    :: dt = 0.001_rk
    integer(ik),         save    :: nsteps = 100
   
   
   
    ! MATRIX SOLVER
    !--------------------------------------------------
    character(len=100),  save    :: msolver = 'direct'
    real(rk),            save    :: tol = 1.e-8
   
   
    ! IO
    !--------------------------------------------------
    integer(ik),         save    :: output_res     = 5
    integer(ik),         save    :: nwrite         = 100
    character(len=100),  save    :: tecplot_prefix = 'tec'
    character(len=100),  save    :: hdf_out        = 'solution.h5'
     
    
       
    ! BOUNDARY CONDITIONS
    !--------------------------------------------------
    integer(ik),         save    :: bc_ximin(MAXBLOCKS),   bc_ximax(MAXBLOCKS)
    integer(ik),         save    :: bc_etamin(MAXBLOCKS),  bc_etamax(MAXBLOCKS)
    integer(ik),         save    :: bc_zetamin(MAXBLOCKS), bc_zetamax(MAXBLOCKS)
    real(rk),            save    :: bcpar1(6), bcpar2(6), bcpar3(6), bcpar4(6)



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



        namelist /polynomial/               basis

        namelist /grid/                     gridfile, gridtype

        namelist /solution/                 solution_order, &
                                            ntime_instances, &
                                            tecplot_prefix,  &
                                            hdf_out 

        namelist /quadrature/               gq_rule

        namelist /equation_set/             eqnset

        namelist /solution_advancement/     temporal_scheme, dt, nsteps

        namelist /matrix_solver/            msolver, tol

!        namelist /grid/     gridfile,   gridtype,   &
!                            bc_ximin,   bc_ximax,   &
!                            bc_etamin,  bc_etamax,  &
!                            bc_zetamin, bc_zetamax, &
!                            bcpar1, bcpar2, bcpar3, bcpar4

        namelist /io/       nwrite, output_res

        inquire(file='chidg.nml', exist=file_exists)
        if (.not. file_exists) call signal(FATAL, "read_input: 'chidg.nml' input file was not found")

        open(unit=7,form='formatted',file="chidg.nml")
        read(7,nml=grid)
        read(7,nml=solution)
        read(7,nml=polynomial)
        read(7,nml=quadrature)
        read(7,nml=equation_set)
        read(7,nml=solution_advancement)
        read(7,nml=matrix_solver)
        read(7,nml=io)



        ! Compute number of terms in polynomial expansions
        nterms_sol1d = (solution_order + 1)
        nterms_s = nterms_sol1d * nterms_sol1d * nterms_sol1d


    end subroutine read_input




end module mod_io
