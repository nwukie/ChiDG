!>  Chimera-based, discontinuous Galerkin equation solver
!!
!!  This program is designed to solve partial differential equations,
!!  and systems of partial differential equations, using the discontinuous
!!  Galerkin method for spatial discretization using Chimera, overset grids to
!!  represent the simulation domain.
!!
!!  @author Nathan A. Wukie
!!  @date   1/31/2016
!!
!!
!---------------------------------------------------------------------------------------------
program driver
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use mod_constants,              only: CHIMERA, NFACES
    use type_chidg,                 only: chidg_t
    use mod_grid_operators,         only: initialize_variable
    use type_function,              only: function_t
    use mod_function,               only: create_function
    use mod_io

    ! Actions
    use mod_chidg_edit,         only: chidg_edit
    use mod_chidg_convert,      only: chidg_convert
    use mod_chidg_interpolate,  only: chidg_interpolate
    use mod_chidg_post,         only: chidg_post
!    use mod_kirchoffs,          only: kirchoff

    ! MPI
    use mod_chidg_mpi,          only: IRANK, NRANK
    use mpi_f08,                only: MPI_COMM_WORLD

    
    !
    ! Variable declarations
    !
    implicit none
    type(chidg_t)                               :: chidg



    integer                                     :: ierr, narg
    integer(ik)                                 :: elems(6), iproc, idom_search, ielem_search, idom, ielem, iface, idonor, nelem_search, ChiID, donor_domain, donor_element
    integer(ik)                                 :: recv_comm, recv_element, recv_domain, donor_proc
    character(len=1024)                         :: chidg_action, filename, file_a, file_b
    class(function_t),              allocatable :: constant, monopole, fcn, polynomial
    logical                                     :: chimera_face






    !
    ! Check for command-line arguments
    !
    narg = command_argument_count()


    !
    ! Execute ChiDG calculation
    !
    if ( narg == 0 ) then



        !
        ! Initialize ChiDG environment
        !
        call chidg%init('mpi')
        call chidg%init('env',MPI_COMM_WORLD)
        call chidg%init('io')



        !
        ! Read grid and boundary condition data
        !
        call chidg%read_grid(gridfile, spacedim)
        call chidg%read_boundaryconditions(gridfile)



        !
        ! Set ChiDG components
        !
        call chidg%set('time_scheme',      time_scheme,      toptions)
        call chidg%set('nonlinear_solver', nonlinear_solver, noptions)
        call chidg%set('linear_solver',    linear_solver,    loptions)
        call chidg%set('preconditioner',   preconditioner)


        !
        ! Initialize solution data storage
        !
        call chidg%initialize_solution_domains(nterms_s)
        call chidg%init('communication')
        call chidg%init('chimera')
        call chidg%initialize_solution_solver()



        !
        ! Initialize solution
        !
        if (solutionfile_in == 'none') then

!            !
!            ! Set initial solution
!            !
!            call create_function(fcn,'gaussian')
!            call fcn%set_option('b_x',0._rk)
!            call fcn%set_option('b_y',1.5_rk)
!            call fcn%set_option('b_z',1.5_rk)
!            call fcn%set_option('c',1.0_rk)
!            call initialize_variable(chidg%data,1,fcn)


!            call create_function(constant,'constant')
!            call create_function(polynomial,'polynomial')
!
!            ! rho
!            call constant%set_option('val',0.01_rk)
!
!            call polynomial%set_option('f',3.5_rk)
!            call initialize_variable(chidg%data,1,polynomial)


            call create_function(constant,'constant')

            ! rho
            call constant%set_option('val',1.19_rk)
            call initialize_variable(chidg%data,1,constant)

            ! rho_u
            call constant%set_option('val',10._rk)
            call initialize_variable(chidg%data,2,constant)

            ! rho_v
            call constant%set_option('val',0._rk)
            call initialize_variable(chidg%data,3,constant)

            ! rho_w
            call constant%set_option('val',0._rk)
            call initialize_variable(chidg%data,4,constant)

            ! rho_E
            call constant%set_option('val',250000._rk)
            call initialize_variable(chidg%data,5,constant)



        else

            !
            ! TODO: put in check that solutionfile actually contains solution
            !
            call chidg%read_solution(solutionfile_in)

        end if

        


        !
        ! Wrap-up initialization activities
        !
        call chidg%init('finalize')

        !
        ! Write initial solution
        !
        if (initial_write) call chidg%write_solution(solutionfile_out)



        !
        ! Run ChiDG simulation
        !
        call chidg%report('before')
        call chidg%run()
        call chidg%report('after')





        !
        ! Write final solution
        !
        if (final_write) call chidg%write_solution(solutionfile_out)



        !
        ! Close ChiDG
        !
        call chidg%close('core')
        call chidg%close('mpi')









    !
    ! ChiDG tool execution. 2 arguments.
    !
    else if ( narg == 2 ) then


        call get_command_argument(1,chidg_action)
        call get_command_argument(2,filename)
        chidg_action = trim(chidg_action)
        filename = trim(filename)
        

        !
        ! Initialize ChiDG environment
        !
        call chidg%init('env')


        !
        ! Select ChiDG action
        !
        if ( trim(chidg_action) == 'edit' ) then
            call chidg_edit(trim(filename))

        else if ( trim(chidg_action) == 'convert' ) then
            call chidg_convert(trim(filename))

        else if ( trim(chidg_action) == 'post' ) then
            call chidg_post(trim(filename))

!        else if ( trim(chidg_action) == 'kirchoff' ) then
!            call kirchoff(filename)


        else
            call chidg_signal(FATAL,"chidg: unrecognized action '"//trim(chidg_action)//"'. Valid options are: 'edit', 'convert'")

        end if


        !
        ! Close ChiDG interface
        !
        call chidg%close('core')

    !
    ! ChiDG tool execution. 3 arguments.
    !
    else if ( narg == 3 ) then


        call get_command_argument(1,chidg_action)
        call get_command_argument(2,file_a)
        call get_command_argument(3,file_b)
        


        !
        ! Initialize ChiDG environment
        !
        call chidg%init('env')



        if ( trim(chidg_action) == 'interpolate' ) then
            call chidg_interpolate(trim(file_a), trim(file_b))

        else
            call chidg_signal(FATAL,"chidg: unrecognized action '"//trim(chidg_action)//"'. Valid options are: 'edit', 'convert'")

        end if



        !
        ! Close ChiDG interface
        !
        call chidg%close('core')




    else
        call chidg_signal(FATAL,"chidg: invalid number of arguments. Expecting (0) arguments: 'chidg'. or (2) arguments: 'chidg action file'.")
    end if












end program driver
