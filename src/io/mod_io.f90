!>  Namelist specifications for IO.
!!
!!  @author Nathan A. Wukie
!!  @date   2/3/2016
!!
!!
!----------------------------------------------------------------------------------------------
module mod_io
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: ZERO
    use type_dict,      only: dict_t
    implicit none

    ! Namelist read flag
    logical                         :: namelist_read_error = .false.

    ! Files
    character(len=100),     save    :: gridfile         = 'none'
    character(len=100),     save    :: solutionfile_in  = 'none'
    character(len=100),     save    :: solutionfile_out = 'none'

    ! Space
    character(len=100),     save    :: basis            = 'legendre'
    integer(ik),            save    :: solution_order   = 1

    ! Quadrature
    integer(ik),            save    :: gq_rule          = 2          ! 1: Collocation, 2: Over-integration
   
    ! Time
    character(len=100),     save    :: time_integrator  = 'steady'
    real(rk),               save    :: dt               = 0.001_rk
    integer(ik),            save    :: time_steps       = 100
    integer(ik),            save    :: ntime_instances  = 1         !TODO: probably we should ask for the entire time of the analysis rather than dt and how many steps
    real(rk),               save    :: frequencies(100) = ZERO
   
    ! Nonlinear solver parameters
    !   ntol    : absolute convergence tolerance
    !   nrtol   : relative convergence tolerance
    !   nnmax   : max number of nonlinear iterations
    !   search  : line-search algorithm ('Backtrack','none')
    !   ptc     : pseudo-transient continuation
    !   smooth  : residual smoothing
    !   cfl0    : initial cfl scaling for pseudo-transient continuation 
    !   cflmax  : max cfl during pseudo-transient continuation. (< 0) = unlimited
    character(len=100),     save    :: nonlinear_solver  = 'newton'
    real(rk),               save    :: ntol              = 1.e-8
    real(rk),               save    :: nrtol             = 1.e-8
    integer(ik),            save    :: nnmax             = -1
    character(len=100),     save    :: search            = 'Backtrack'
    logical,                save    :: ptc               = .true.
    logical,                save    :: smooth            = .true.
    real(rk),               save    :: cfl0              = 1._rk
    real(rk),               save    :: cflmax            = -1._rk
    type(dict_t),           save    :: noptions
    
    ! Linear solver parameters
    !   preconditioner          : 'Identity','Jacobi','ILU0','RASILU0'
    !   ltol                    : absolute convergence tolerance
    !   lrtol                   : relative convergence tolerance
    !   lnmax                   : max number of linear iterations
    !   nkrylov                 : number of krylov vectors
    !   orthogonalization       : 'CGS','MGS'
    !   inner_fgmres            : inner fgmres iterative preconditioning
    !   inner_ltol              : absolute convergence tolerance for inner iteration
    !   inner_lrtol             : relative convergence tolerance for inner iteration
    !   inner_lnmax             : maximum number of inner iterations
    !   inner_nkrylov           : number of inner krylov vectors
    !   inner_orthogonalization : 'CGS','MGS'
    !   inner_silence           : -10(no report of inner iteration), 0(report inner iteration)
    character(len=100),     save    :: linear_solver    = 'fgmres'
    character(len=100),     save    :: preconditioner   = 'identity'
    real(rk),               save    :: ltol             = 1.e-8_rk
    real(rk),               save    :: lrtol            = 1.e-8_rk
    integer(ik),            save    :: lnmax            = -1
    integer(ik),            save    :: nkrylov          = 2000
    character(len=3),       save    :: orthogonalization = 'CGS'
    logical,                save    :: inner_fgmres     = .true.
    real(rk),               save    :: inner_ltol       = 1.e-1_rk
    real(rk),               save    :: inner_lrtol      = 1.e-1_rk
    integer(ik),            save    :: inner_lnmax      = 100
    integer(ik),            save    :: inner_nkrylov    = 100
    character(len=3),       save    :: inner_orthogonalization = 'CGS'
    integer(ik),            save    :: inner_silence    = -10
    type(dict_t),           save    :: loptions
   
    ! Output
    !   initial_write   : write initial solution
    !   final_write     : write final solution
    !   nwrite          : write solution every 'nwrite' steps of time-integrator
    !   verbosity       : 1-5, controls amount of information to screen/log file
    logical,                save    :: initial_write    = .false.
    logical,                save    :: final_write      = .true.
    integer(ik),            save    :: nwrite           = 100
    integer(ik),            save    :: verbosity        = 2
     

    ! Initial fields
    !
    !   Specify constant initial solution. 
    !
    !   READ ONLY IF: solutionfile_in='none'
    !
    !   Example for Euler fluid problem(5 equations):
    !   ---------------------------------------------
    !   &initial
    !       initial_fields = 1.19, 150., 0., 0., 250000.
    !   /
    real(rk),               save    :: initial_fields(100) = ZERO


    ! Namelist Groups
    namelist /files/                    gridfile,           &
                                        solutionfile_in,    &
                                        solutionfile_out

    namelist /space/                    basis,              &
                                        solution_order

    namelist /quadrature/               gq_rule


    namelist /time/                     time_integrator,    &
                                        dt,                 &
                                        time_steps,         &
                                        ntime_instances,    &
                                        frequencies


    namelist /nonlinear_solve/          nonlinear_solver,   &
                                        ntol,               &
                                        nrtol,              &
                                        nnmax,              &
                                        cfl0,               &
                                        cflmax,             &
                                        search,             &
                                        ptc,                &
                                        smooth

    namelist /linear_solve/             linear_solver,      &
                                        ltol,               &
                                        lrtol,              &
                                        lnmax,              &
                                        preconditioner,     &
                                        nkrylov,            &
                                        orthogonalization,  &
                                        inner_fgmres,       &
                                        inner_ltol,         &
                                        inner_lrtol,        &
                                        inner_lnmax,        &
                                        inner_nkrylov,      &
                                        inner_silence,      &
                                        inner_orthogonalization


    namelist /io/                       nwrite,             &
                                        initial_write,      &
                                        final_write,        &
                                        verbosity

    namelist /initial/                  initial_fields



contains



    !>  Initialize solver dictionaries with default values.
    !!
    !!  Called from chidg%start_up('core'). In case 'read_input' isn't called
    !!  (for example, in a test) the dictionaries are still initialized so the 
    !!  initialization of linear and nonlinear solvers are able to access the 
    !!  default parameters.
    !!
    !!  @author Nathan A. Wukie (AFRL)
    !!  @date   12/28/2018
    !!
    !-----------------------------------------------------------------------------
    subroutine initialize_input_dictionaries(noptions_in,loptions_in)
        type(dict_t)    :: noptions_in
        type(dict_t)    :: loptions_in

        ! Set nonlinear solver options
        call noptions_in%set('tol',ntol)
        call noptions_in%set('rtol',nrtol)
        call noptions_in%set('cfl0',cfl0)
        call noptions_in%set('cflmax',cflmax)
        call noptions_in%set('nmax',nnmax)
        call noptions_in%set('search',search)
        call noptions_in%set('ptc',ptc)
        call noptions_in%set('smooth',smooth)

        ! Set linear solver options
        call loptions_in%set('tol',ltol)
        call loptions_in%set('rtol',lrtol)
        call loptions_in%set('nmax',lnmax)
        call loptions_in%set('nkrylov',nkrylov)
        call loptions_in%set('orthogonalization',orthogonalization)

        call loptions_in%set('inner_fgmres',inner_fgmres)
        call loptions_in%set('inner_tol',inner_ltol)
        call loptions_in%set('inner_rtol',inner_lrtol)
        call loptions_in%set('inner_nmax',inner_lnmax)
        call loptions_in%set('inner_nkrylov',inner_nkrylov)
        call loptions_in%set('inner_silence',inner_silence)
        call loptions_in%set('inner_orthogonalization',inner_orthogonalization)


    end subroutine initialize_input_dictionaries
    !*****************************************************************************




    !> Read input file
    !!      - Read chidg.nml namelist formatted input file
    !!
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/3/2016
    !!
    !-------------------------------------------------------------------------------------------
    subroutine read_input()
        logical :: file_exists
        integer :: file_unit, msg

        ! Check that input file exists
        inquire(file='chidg.nml', exist=file_exists)
        if (.not. file_exists) call chidg_signal(FATAL, "read_input: 'chidg.nml' input file was not found")


        ! Read namelist input for parameter initialization
!        open(newunit=file_unit,form='formatted',file="chidg.nml")
!        read(file_unit,nml=files,           iostat=msg)
!        if (msg /= 0) call chidg_signal(FATAL,"read_input: error reading 'files' namelist from 'chidg.nml'.")
!        read(file_unit,nml=space,           iostat=msg)
!        if (msg /= 0) call chidg_signal(FATAL,"read_input: error reading 'space' namelist from 'chidg.nml'.")
!        read(file_unit,nml=quadrature,      iostat=msg)
!        if (msg /= 0) call chidg_signal(FATAL,"read_input: error reading 'quadrature' namelist from 'chidg.nml'.")
!        read(file_unit,nml=time,            iostat=msg)
!        if (msg /= 0) call chidg_signal(FATAL,"read_input: error reading 'time' namelist from 'chidg.nml'.")
!        read(file_unit,nml=nonlinear_solve, iostat=msg)
!        if (msg /= 0) call chidg_signal(FATAL,"read_input: error reading 'nonlinear_solve' namelist from 'chidg.nml'.")
!        read(file_unit,nml=linear_solve,    iostat=msg)
!        if (msg /= 0) call chidg_signal(FATAL,"read_input: error reading 'linear_solve' namelist from 'chidg.nml'.")
!        read(file_unit,nml=io,              iostat=msg)
!        if (msg /= 0) call chidg_signal(FATAL,"read_input: error reading 'io' namelist from 'chidg.nml'.")
!        read(file_unit,nml=initial,         iostat=msg)
!        if (msg /= 0) call chidg_signal(FATAL,"read_input: error reading 'initial' namelist from 'chidg.nml'.")
!        close(file_unit)

        open(newunit=file_unit,form='formatted',file="chidg.nml")
        read(file_unit,nml=files,           iostat=msg)
        if (msg/=0) call handle_namelist_error(file_unit,msg)
        read(file_unit,nml=space,           iostat=msg)
        if (msg/=0) call handle_namelist_error(file_unit,msg)
        read(file_unit,nml=quadrature,      iostat=msg)
        if (msg/=0) call handle_namelist_error(file_unit,msg)
        read(file_unit,nml=time,            iostat=msg)
        if (msg/=0) call handle_namelist_error(file_unit,msg)
        read(file_unit,nml=nonlinear_solve, iostat=msg)
        if (msg/=0) call handle_namelist_error(file_unit,msg)
        read(file_unit,nml=linear_solve,    iostat=msg)
        if (msg/=0) call handle_namelist_error(file_unit,msg)
        read(file_unit,nml=io,              iostat=msg)
        if (msg/=0) call handle_namelist_error(file_unit,msg)
        read(file_unit,nml=initial,         iostat=msg)
        if (msg/=0) call handle_namelist_error(file_unit,msg)
        close(file_unit)

        ! Exit if error detected in namelist read
        if (namelist_read_error) call chidg_abort()





        ! Set nonlinear solver options
        call noptions%set('tol',ntol)
        call noptions%set('rtol',nrtol)
        call noptions%set('nmax',nnmax)
        call noptions%set('cfl0',cfl0)
        call noptions%set('cflmax',cflmax)

        ! Set linear solver options
        call loptions%set('tol',ltol)
        call loptions%set('rtol',lrtol)
        call loptions%set('nmax',lnmax)
        call loptions%set('nkrylov',nkrylov)
        call loptions%set('orthogonalization',orthogonalization)

        call loptions%set('inner_fgmres',inner_fgmres)
        call loptions%set('inner_tol',inner_ltol)
        call loptions%set('inner_rtol',inner_lrtol)
        call loptions%set('inner_nmax',inner_lnmax)
        call loptions%set('inner_nkrylov',inner_nkrylov)
        call loptions%set('inner_silence',inner_silence)
        call loptions%set('inner_orthogonalization',inner_orthogonalization)


    end subroutine read_input
    !****************************************************************************************




    !>  Handle namelist read error.
    !!
    !!  @author Nathan A. Wukie(AFRL)
    !!  @date   1/3/2019
    !!
    !!  @attribution Jacob Williams
    !!  @license     See label: JWILLIAMS in LICENSE file
    !!  @purpose     Approach for namelist error checking
    !!
    !---------------------------------------------------------------------------------------
    subroutine handle_namelist_error(file_unit,msg)
        use iso_fortran_env,    only: error_unit
        integer,        intent(in)  :: file_unit 
        integer,        intent(in)  :: msg

        character(len=1000) :: line

        backspace(file_unit)
        read(file_unit,fmt='(A)') line
        write(error_unit,'(A)') 'Invalid line in namelist: '//trim(line) 
        namelist_read_error = .true.

    end subroutine handle_namelist_error
    !***************************************************************************************





    !>  Write a file containing the default namelist entries.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   5/17/2017
    !!
    !----------------------------------------------------------------------------------------
    subroutine write_namelist()
        integer :: file_unit

        ! Write default namelist input for parameter initialization
        open(newunit=file_unit,form='formatted',file="chidg.nml")
        write(unit=file_unit, nml=files)
        write(unit=file_unit, nml=space)
        write(unit=file_unit, nml=quadrature)
        write(unit=file_unit, nml=time)
        write(unit=file_unit, nml=nonlinear_solve)
        write(unit=file_unit, nml=linear_solve)
        write(unit=file_unit, nml=io)
        write(unit=file_unit, nml=initial)
        close(unit=file_unit)

    end subroutine write_namelist
    !****************************************************************************************



end module mod_io
