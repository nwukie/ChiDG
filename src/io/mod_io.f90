!>  Namelist specifications for IO.
!!
!!  @author Nathan A. Wukie
!!  @date   2/3/2016
!!
!!
!----------------------------------------------------------------------------------------------------------
module mod_io
#include <messenger.h>
    !!          Module Includes             !!
    use mod_kinds,     only: rk,ik
    use mod_constants, only: MAXBLOCKS, TWO_DIM, THREE_DIM

    !!          Variable Declarations       !!
    implicit none


    ! FILES
    !--------------------------------------------------
    character(len=100),  save    :: gridfile
    character(len=100),  save    :: gridtype
    character(len=100),  save    :: tecplot_prefix   = 'tec'
    character(len=100),  save    :: hdf_out          = 'solution.h5'

    character(len=100),  save    :: solutionfile_in  = 'none'
    character(len=100),  save    :: solutionfile_out = 'none'




    ! SPACE
    !--------------------------------------------------
    character(len=100),  save    :: basis            = 'legendre'
    integer(ik),         save    :: solution_order   = 1
    integer(ik),         save    :: spacedim         = 3

 


    
    ! QUADRATURE
    !--------------------------------------------------
    integer(ik),         save    :: gq_rule          = 2          !> 1: Collocation, 2: Over-integration
   


   
    
    ! EQUATION SET 
    !--------------------------------------------------
    character(len=100),  save    :: eqnset           = 'scalar'
  
  



    
    ! TIME
    !--------------------------------------------------
    character(len=100),  save    :: time_scheme      = 'steady'
    real(rk),            save    :: dt               = 0.001_rk
    integer(ik),         save    :: time_steps       = 100
    real(rk),            save    :: ttol             = 1.e-8
    integer(ik),         save    :: ntime_instances  = 1
   
   

    ! NONLINEAR SOLVER
    !-------------------------------------------------
    character(len=100),  save    :: nonlinear_solver = 'newton'
    integer(ik),         save    :: nonlinear_steps  = 100
    real(rk),            save    :: ntol             = 1.e-8
    real(rk),            save    :: cfl0             = 1._rk
    

   
    ! LINEAR SOLVER
    !--------------------------------------------------
    character(len=100),  save    :: linear_solver    = 'fgmres'
    real(rk),            save    :: ltol             = 1.e-8



    ! PRECONDITIONER
    !--------------------------------------------------
    character(len=100),  save    :: preconditioner   = 'identity'
   
   



    ! IO
    !--------------------------------------------------
    integer(ik),         save    :: nwrite           = 100
    logical,             save    :: initial_write    = .true.
    logical,             save    :: final_write      = .true.
    integer(ik),         save    :: output_res       = 10
     
    



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
    !!  @date   2/3/2016
    !!
    !----------------------------------------------------------------------------------------------------------
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
                                            solution_order,        &
                                            spacedim

        namelist /quadrature/               gq_rule

        namelist /equation_set/             eqnset

        namelist /time/                     time_scheme,           &
                                            dt,                    &
                                            time_steps,            &
                                            ntime_instances,       &
                                            ttol


        namelist /nonlinear_solve/          nonlinear_solver,      &
                                            nonlinear_steps,       &
                                            cfl0,                  &
                                            ntol

        namelist /linear_solve/             linear_solver,         &
                                            ltol,                  &
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
        read(7,nml=nonlinear_solve)
        read(7,nml=linear_solve)
        read(7,nml=io)



        !
        ! Compute number of terms in polynomial expansions
        !
        nterms_sol1d = (solution_order)

        if ( spacedim == THREE_DIM ) then
            nterms_s = nterms_sol1d * nterms_sol1d * nterms_sol1d
        else if ( spacedim == TWO_DIM ) then
            nterms_s = nterms_sol1d * nterms_sol1d
        else
            call chidg_signal(FATAL,"mod_io: Invalid spacedim")
        end if



    end subroutine read_input
    !*******************************************************************************************************




end module mod_io
