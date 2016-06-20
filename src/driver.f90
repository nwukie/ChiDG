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
    use mod_constants,              only: ONE, TWO, ZERO
    use type_chidg,                 only: chidg_t
    use mod_grid_operators,         only: initialize_variable
    use type_dict,                  only: dict_t
    use type_function,              only: function_t
    use mod_function,               only: create_function
    use type_domain_connectivity,   only: domain_connectivity_t
    use type_partition,             only: partition_t
    use mod_partitioners,           only: partition_connectivity, send_partitions, recv_partition
    use mod_hdfio,                  only: read_connectivity_hdf
    use mod_io

    ! Actions
    use mod_chidg_edit,         only: chidg_edit
    use mod_chidg_convert,      only: chidg_convert
    use mod_chidg_interpolate,  only: chidg_interpolate
    use mod_chidg_post,         only: chidg_post
    use mod_kirchoffs,          only: kirchoff

    ! MPI
    use mod_chidg_mpi,          only: IRANK, NRANK, GLOBAL_MASTER
    use mpi_f08,                only: MPI_Barrier, MPI_COMM_WORLD

    
    !
    ! Variable declarations
    !
    implicit none
!    include 'mpif.h'
    type(chidg_t)                       :: chidg
    type(partition_t)                   :: partition
    type(dict_t)                        :: toptions, noptions, loptions

    type(domain_connectivity_t),    allocatable :: connectivities(:)
    type(partition_t),              allocatable :: partitions(:)

    class(function_t),  allocatable     :: constant, monopole

    integer                             :: ierr, narg, iread
    character(len=1024)                 :: chidg_action, filename, file_a, file_b






    !
    ! Check for command-line arguments
    !
    narg = command_argument_count()


    !
    ! Execute ChiDG calculation
    !
    if ( narg == 0 ) then


!        !
!        ! Initialize MPI
!        !
!        call MPI_Init(ierr)
!        if (ierr /= 0) call chidg_signal(FATAL,"MPI_Init")
!        call MPI_Comm_Size(MPI_COMM_WORLD,NRANK,ierr)
!        if (ierr /= 0) call chidg_signal(FATAL,"MPI_Comm_Size")
!        call MPI_Comm_Rank(MPI_COMM_WORLD,IRANK,ierr)
!        if (ierr /= 0) call chidg_signal(FATAL,"MPI_Comm_Rank")



        !
        ! Initialize ChiDG environment
        !
        call chidg%init('mpi')
        call chidg%init('env')
        call chidg%init('io')



        !
        ! Read connectivity, partition, distribute
        !
        if ( irank == GLOBAL_MASTER ) call write_line("Partitioning grid")
        if ( irank == GLOBAL_MASTER ) then
            call read_connectivity_hdf(gridfile,connectivities)

            call partition_connectivity(connectivities, partitions)

            call send_partitions(partitions)
        end if



        !
        ! Receive partition from GLOBAL_MASTER
        !
        call recv_partition(partition)




        !
        ! Read grid data from file. One partition at a time to avoid the file being accessed simultaneously.
        !
        if ( irank == GLOBAL_MASTER ) call write_line("Reading grid")
        do iread = 0,NRANK-1
            if ( iread == IRANK ) then
                print*, "Rank ", IRANK, ": reading grid"
                call chidg%read_grid(gridfile, spacedim, partition)
                print*, "Rank ", IRANK, ": done reading grid"
                !print*, IRANK, " reading boundary conditions"
                !call chidg%read_boundaryconditions(gridfile, partition)
            end if

            print*, "Rank ", IRANK, ": waiting"
            call MPI_Barrier(MPI_COMM_WORLD,ierr)  ! sync to prevent simultaneous file access
        end do






        ! Set time-scheme options
        call toptions%set('dt',dt)
        call toptions%set('nsteps',time_steps)
        call toptions%set('nwrite',nwrite)


        ! Set nonlinear solver options
        call noptions%set('tol',ntol)
        call noptions%set('cfl0',cfl0)
        call noptions%set('nsteps',nonlinear_steps)

        ! Set linear solver options
        call loptions%set('tol',ltol)



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
        print*, "Rank ", IRANK, ": Initializing grid solution data"
        call chidg%initialize_solution_domains(nterms_s)
        print*, "Rank ", IRANK, ": Initializing communication"
        call chidg%init('communication')
!        call chidg%init('chimera')
!        call chidg%initialize_solution_solver()
!
!
!
!        !
!        ! Initialize solution
!        !
!        if (solutionfile_in == 'none') then
!            call create_function(constant,'constant')
!
!
!            ! rho
!            call constant%set_option('val',1.20_rk)
!            call initialize_variable(chidg%data,1,constant)
!
!            ! rho_u
!            call constant%set_option('val',50._rk)
!            call initialize_variable(chidg%data,2,constant)
!
!            ! rho_v
!            call constant%set_option('val',0._rk)
!            call initialize_variable(chidg%data,3,constant)
!
!            ! rho_w
!            call constant%set_option('val',0._rk)
!            call initialize_variable(chidg%data,4,constant)
!
!            ! rho_E
!            call constant%set_option('val',230000._rk)
!            call initialize_variable(chidg%data,5,constant)
!
!
!!            ! rho
!!            call constant%set_option('val',0._rk)
!!            call initialize_variable(chidg%data,6,constant)
!!
!!            ! rho_u
!!            call constant%set_option('val',0._rk)
!!            call initialize_variable(chidg%data,7,constant)
!!
!!            ! rho_v
!!            call constant%set_option('val',0._rk)
!!            call initialize_variable(chidg%data,8,constant)
!!
!!            ! rho_w
!!            call constant%set_option('val',0._rk)
!!            call initialize_variable(chidg%data,9,constant)
!!
!!            ! rho_E
!!            call constant%set_option('val',0._rk)
!!            call initialize_variable(chidg%data,10,constant)
!!
!
!        else
!
!            !
!            ! TODO: put in check that solutionfile actually contains solution
!            !
!            call chidg%read_solution(solutionfile_in)
!!            do iread = 1,nrank
!!                if ( iread == irank ) then
!!                    call chidg%read_solution(solutionfile_in,partition)
!!                end if
!!            end do
!
!        end if
!
!        
!
!
!        !
!        ! Wrap-up initialization activities
!        !
!        call chidg%init('finalize')
!
!        !
!        ! Write initial solution
!        !
!        if (initial_write) call chidg%write_solution(solutionfile_out)
!
!
!
!
!
!        !
!        ! Run ChiDG simulation
!        !
!        call chidg%run()
!
!
!
!
!
!        !
!        ! Write final solution
!        !
!        if (final_write) call chidg%write_solution(solutionfile_out)
!
!        !
!        ! Reporting
!        !
!        call chidg%report()
!
!!        !
!!        ! Close MPI
!!        !
!!        call MPI_Finalize(ierr)





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

        else if ( trim(chidg_action) == 'kirchoff' ) then
            call kirchoff(filename)


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
