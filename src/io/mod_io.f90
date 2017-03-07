!>  Namelist specifications for IO.
!!
!!  @author Nathan A. Wukie
!!  @date   2/3/2016
!!
!!
!----------------------------------------------------------------------------------------------------------
module mod_io
#include <messenger.h>
    use mod_kinds,      only: rk,ik
    use mod_constants,  only: MAXBLOCKS, TWO_DIM, THREE_DIM, ZERO
    use type_dict,      only: dict_t
    implicit none


    ! Files
    character(len=100),     save    :: gridfile
    character(len=100),     save    :: gridtype
    character(len=100),     save    :: solutionfile_in  = 'none'
    character(len=100),     save    :: solutionfile_out = 'none'

    ! Space
    character(len=100),     save    :: basis            = 'legendre'
    integer(ik),            save    :: solution_order   = 1
    integer(ik),            save    :: spacedim         = 3

    ! Quadrature
    integer(ik),            save    :: gq_rule          = 2          !> 1: Collocation, 2: Over-integration
   
    ! Time
    character(len=100),     save    :: time_integrator  = 'steady'
    real(rk),               save    :: dt               = 0.001_rk
    integer(ik),            save    :: time_steps       = 100
    real(rk),               save    :: ttol             = 1.e-8     !Apparently not used
    integer(ik),            save    :: ntime_instances  = 1         !TODO: probably we should ask for the entire time of the analysis rather than dt and how many steps
    real(rk),               save    :: frequencies(100) = ZERO
!    type(dict_t),           save    :: toptions
   
    ! Nonlinear solver
    character(len=100),     save    :: nonlinear_solver = 'newton'
    integer(ik),            save    :: nonlinear_steps  = 100
    real(rk),               save    :: ntol             = 1.e-8
    real(rk),               save    :: cfl0             = 1._rk
    type(dict_t),           save    :: noptions
    
    ! Linear solver
    character(len=100),     save    :: linear_solver    = 'fgmres'
    real(rk),               save    :: ltol             = 1.e-8
    character(len=100),     save    :: preconditioner   = 'identity'
    type(dict_t),           save    :: loptions
   
    ! io
    integer(ik),            save    :: nwrite           = 100
    logical,                save    :: initial_write    = .false.
    logical,                save    :: final_write      = .true.
!    integer(ik),         save    :: output_res       = 10
     



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
                                            solutionfile_in,       &
                                            solutionfile_out

        namelist /space/                    basis,                 &
                                            solution_order,        &
                                            spacedim

        namelist /quadrature/               gq_rule


        namelist /time/                     time_integrator,       &
                                            dt,                    &
                                            time_steps,            &
                                            ttol,                  &
                                            ntime_instances,       &
                                            frequencies


        namelist /nonlinear_solve/          nonlinear_solver,      &
                                            nonlinear_steps,       &
                                            cfl0,                  &
                                            ntol

        namelist /linear_solve/             linear_solver,         &
                                            ltol,                  &
                                            preconditioner


        namelist /io/                       nwrite,                &
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
        read(7,nml=time)
        read(7,nml=nonlinear_solve)
        read(7,nml=linear_solve)
        read(7,nml=io)




        !
        ! Initialize options dictionaries
        !


        ! Set time-integrator options
!        call toptions%set('dt',dt)                 !
!        call toptions%set('nsteps',time_steps)     ! TO BE REMOVED
!        call toptions%set('nwrite',nwrite)         !

        ! Set nonlinear solver options
        call noptions%set('tol',ntol)
        call noptions%set('cfl0',cfl0)
        call noptions%set('nsteps',nonlinear_steps)

        ! Set linear solver options
        call loptions%set('tol',ltol)



    end subroutine read_input
    !*******************************************************************************************************




end module mod_io
