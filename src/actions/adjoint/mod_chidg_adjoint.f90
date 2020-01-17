
!!  adjoint variables for a given objective function previously defined in the non-linear(CFD) 
!!  solver.
!!  It reads in the solution file written by the non-linear computation and solves the adjoint
!!  linear system of equations.
!!  It uses the linear solver and preconditioner specified in the namelist (chidg.nml). 
!!  The adjoint variables are written in the same solution file of the non-linear solver
!!
!!  @author Matteo Ugolotti
!!  @date 9/12/2017
!!
!!
!!  Usage: chidg adjoint
!!  
!
!----------------------------------------------------------------------------------------------
module mod_chidg_adjoint
#include <messenger.h>
    use mod_kinds,                  only: rk, ik
    use type_chidg,                 only: chidg_t
    use mod_chidg_mpi,              only: GLOBAL_MASTER, ChiDG_COMM, IRANK, NRANK
    use type_svector,               only: svector_t
    use mod_string,                 only: string_t, get_file_prefix, check_file_extension
    use type_file_properties,       only: file_properties_t
    use mod_hdf_utilities,          only: get_properties_hdf
    use mod_file_utilities,         only: delete_file
    use mod_io
    implicit none




contains



    !>  Driver for solving adjoint linear system of equations
    !!
    !!  @author Matteo Ugolotti
    !!  @date 9/12/2017
    !!
    !!  Restructured for unsteady adjoint
    !!
    !!  WARNING:
    !!  For unsteady adjoint computation only the unsteady time solution will be considered for
    !!  adjoint variables calculation. Therefore, make sure the initial file for the primal unsteady
    !!  flow solver is named with prefix different from the other step flow files
    !!
    !!  @author Matteo Ugolotti
    !!  @date 8/10/2018
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_adjoint()

        type(chidg_t)                               :: chidg
        integer(ik)                                 :: ierr, ifile_f, istep, iproc
        character(:),                   allocatable :: usr_msg, tmp_flow_file,                  &
                                                       command_flow, prefix, flow_solution_files
        character(len=1024)                         :: flow_solution_file
        character(len=10)                           :: time_string
        type(svector_t)                             :: flow_files
        type(string_t)                              :: str_flow
        type(file_properties_t)                     :: file_props
        logical                                     :: nsteps_mismatch, hdf_file


        ! Initialize ChiDG environment
        call chidg%start_up('mpi')
        call chidg%start_up('namelist')
        call chidg%start_up('core')


        ! Work out wildcards
        prefix = get_file_prefix(solutionfile_in,'.h5')
        flow_solution_files = trim(prefix)//'*'
        !call date_and_time(time=time_string)
        tmp_flow_file = 'chidg_adjoint_flow_files.txt' 
        if (IRANK == GLOBAL_MASTER) then
            command_flow  = 'ls '//trim(flow_solution_files)//' > '//tmp_flow_file
            call system(command_flow) 
        end if
       

        ! Save all files' names in an svector
        do iproc = 0,NRANK-1
            if (iproc == IRANK) then
                open(7,file=tmp_flow_file,action='read')
                do
                    read(7,fmt='(a)', iostat=ierr) flow_solution_file
                    if (ierr /= 0) exit
                    hdf_file = check_file_extension(trim(flow_solution_file),'.h5')
                    if (hdf_file) then
                        call str_flow%set(trim(flow_solution_file))
                        call flow_files%push_back(str_flow)
                    end if
                end do
                close(7)
            end if
            call MPI_Barrier(ChiDG_COMM,ierr)
        end do


        ! Set ChiDG Algorithms, Accuracy
        call chidg%set('Time Integrator' , algorithm=time_integrator                    )
        call chidg%set('Nonlinear Solver', algorithm=nonlinear_solver, options=noptions ) ! Even though not used
        call chidg%set('Linear Solver'   , algorithm=linear_solver,    options=loptions )
        call chidg%set('Preconditioner'  , algorithm=preconditioner                     )
        call chidg%set('Solution Order'  , integer_input=solution_order                 )
        
        call MPI_Barrier(ChiDG_COMM,ierr)

        
        
        ! Theoretically, the number of files should be equal to the number of time steps (nsteps)
        ! defined in the namelist as time_steps
        ! Check that this is true. q_time has been initialized with the number of steps defined 
        ! in the namelist like all the other containers for adjoint and adjointx.
        nsteps_mismatch = ( (flow_files%size() /= time_steps) ) 
        if (nsteps_mismatch) then
            usr_msg = "The number of primal solutions does not match the number of time steps &
                       defined in the namelist. Please, make sure all the files are located in &
                       the current folder or that the time_steps input in the namelist is correct."
            call chidg_signal(FATAL,trim(usr_msg))
        end if



        ! Read grid from solution_file, boundary condition data and functionals
        ! The assumption here is that the primal solver has been run already with functionals
        ! selected or the grid hdf file edited (chidg edit) for adding functionals. 
        call chidg%read_mesh(flow_files%data_(1)%get(),storage='adjoint storage')


        ! Read primal solution/s
        do ifile_f = 1,flow_files%size()
            
            
            ! Get necessary file properties
            do iproc = 0,NRANK-1
                if (iproc == IRANK) then
                    file_props  = get_properties_hdf(flow_files%data_(ifile_f)%get())
                    istep       = file_props%istep
                end if
                call MPI_Barrier(ChiDG_COMM,ierr)
            end do
            
            ! Read in primal solution
            call chidg%read_fields(flow_files%data_(ifile_f)%get(),'primary')
            
            ! Copy primal solution from sdata%q_in to adjoint%q_time(:)
            call chidg%data%sdata%adjoint%process_primal_solution(chidg%data%sdata%q_in,istep)
        
        end do

        ! Run ChiDG Adjoint simulation
        call chidg%reporter('before')
        call chidg%run_adjoint(write_final=final_write)
        call chidg%reporter('after adjoint')


        ! Close ChiDG
        call chidg%shut_down('core')
        call chidg%shut_down('mpi')


        ! Delete tmp files
        if (IRANK == GLOBAL_MASTER) call delete_file(tmp_flow_file)


    end subroutine chidg_adjoint
    !******************************************************************************************




















end module mod_chidg_adjoint
