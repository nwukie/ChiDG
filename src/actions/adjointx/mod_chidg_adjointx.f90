!>  ChiDG AdjointX: this is driver for computing the sensitivities of an objective function
!!  wrt the grid geoemtric variables, that is the grid nodes.
!!  This is an adjoint-based procedure, as such it needs both a primal solution and an adjoint 
!!  solution.
!!  As a ChiDG action it requires one or two arguments: ONE if the file provided contains both 
!!  the primal and adjoint solution, TWO files if primal and adjoint solution are store in two
!!  separate files. Notice that the grid information have to be in the first file provided, 
!!  that is the primal solution.
!!
!!  @author Matteo Ugolotti
!!  @date 7/25/2018
!!
!!
!!  Usage: chidg adjointx
!!  
!
!----------------------------------------------------------------------------------------------
module mod_chidg_adjointx
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



    !>  Driver for computing the grid nodes sensitivities of an objective function
    !!
    !!  @author Matteo Ugolotti
    !!  @date 7/25/2018
    !!
    !!  Required computation
    !!
    !!      DJ/DX = dJ/dX + v^T*dR/dX 
    !!
    !------------------------------------------------------------------------------------------
    subroutine chidg_adjointx()

        type(chidg_t)                               :: chidg
        integer(ik)                                 :: ierr, ifile_f, ifile_a, istep, nfunctionals, &
                                                       iproc
        character(len=:),               allocatable :: usr_msg, tmp_flow_file, tmp_adj_file,    &
                                                       command_flow, command_adj, flow_prefix,  &
                                                       adj_prefix, flow_solution_files,         &
                                                       adjoint_solution_files
        character(len=1024)                         :: flow_solution_file, adjoint_solution_file
        type(svector_t)                             :: flow_files, adjoint_files
        type(string_t)                              :: str_flow, str_adj
        type(file_properties_t)                     :: file_props
        logical                                     :: solutions_mismatch, nsteps_mismatch, hdf_file
        character(len=10)                           :: time_string

        !
        ! Initialize ChiDG environment
        !
        call chidg%start_up('mpi')
        call chidg%start_up('namelist')
        call chidg%start_up('core')


        !
        ! Work out wildcards
        !
        flow_prefix = get_file_prefix(solutionfile_in,'.h5')
        flow_solution_files = trim(flow_prefix)//'*'
        adj_prefix = get_file_prefix(adjointfile_out,'.h5')
        adjoint_solution_files = trim(adj_prefix)//'*'

        !call date_and_time(time=time_string)
        tmp_flow_file = 'chidg_adjointx_flow_files.txt' 
        tmp_adj_file  = 'chidg_adjointx_adj_files.txt' 
        if (IRANK == GLOBAL_MASTER) then
            command_flow  = 'ls '//trim(flow_solution_files)   //' > '//tmp_flow_file
            command_adj   = 'ls '//trim(adjoint_solution_files)//' > '//tmp_adj_file
            call system(command_flow) 
            call system(command_adj) 
        end if
       
        !
        ! Save all files in svector
        !
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
                open(8,file=tmp_adj_file,action='read')
                do
                    read(8,fmt='(a)', iostat=ierr) adjoint_solution_file
                    if (ierr /= 0) exit
                    hdf_file = check_file_extension(trim(adjoint_solution_file),'.h5')
                    if (hdf_file) then
                        call str_adj%set(trim(adjoint_solution_file))
                        call adjoint_files%push_back(str_adj)
                    end if
                end do
                close(8)
            end if
            call MPI_Barrier(ChiDG_COMM,ierr)
        end do
        

        ! Check that the number of primal solutions is equal to number of adjoint solutions
        solutions_mismatch = ( flow_files%size() /= adjoint_files%size() )
        if (solutions_mismatch) then
            call chidg_signal(FATAL,'The number of primal solutions does not match the number of adjoint solutions read in. Make sure that all the files containing the primal and dual solutions are present in the folder. Also check that the file names provided as input argument to ChiDG AdjointX are correct. ')
        end if



        ! Set ChiDG Algorithms, Accuracy
        ! TODO: this might not be necessary
        call chidg%set('Time Integrator' , algorithm=time_integrator                    )
        call chidg%set('Nonlinear Solver', algorithm=nonlinear_solver, options=noptions ) ! not used
        call chidg%set('Linear Solver'   , algorithm=linear_solver,    options=loptions ) ! not used
        call chidg%set('Preconditioner'  , algorithm=preconditioner                     ) ! not used
        call chidg%set('Solution Order'  , integer_input=solution_order                 )
        
        call MPI_Barrier(ChiDG_COMM,ierr)


        ! Theoretically, the number of files should be equal to the number of time steps (nsteps)
        ! defined in the namelist as time_steps
        ! Check that this is true. q_time hs been initialized with the number of steps defined 
        ! in the namelist like all the other containers for adjoint and adjointx.
        nsteps_mismatch = ( (flow_files%size() /= time_steps) .or. (adjoint_files%size() /= time_steps) ) 
        if (nsteps_mismatch) then
            call chidg_signal(FATAL,'The number of primal or adjoint solutions does not match the number of time steps defined in the namelist. Please, make sure all the files are located in the current folder or that the time_steps input in the namelist is correct.')
        end if


        ! Read grid from solution_file, boundary condition data and functionals
        ! The assumption here is that the primal solver has been run already with functionals
        ! selected or edited (chidg edit) for adding functionals. 
        call chidg%read_mesh(flow_files%data_(1)%get(), storage='adjointx storage')

        
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
        
         
        ! Read adjoint solution/s
        do ifile_a = 1,adjoint_files%size()
            
            ! Get necessary file properties
            do iproc = 0,NRANK-1
                if (iproc == IRANK) then
                    file_props   = get_properties_hdf(adjoint_files%data_(ifile_a)%get())
                    istep        = file_props%istep
                    nfunctionals = file_props%nfunctionals
                end if
                call MPI_Barrier(ChiDG_COMM,ierr)
            end do
            
            ! Set adjoint fields and read in adjoint solution
            call chidg%data%set_fields_adjoint(nfunctionals)            
            call chidg%read_fields(adjoint_files%data_(ifile_a)%get(),'adjoint')
            
            ! Copy adjoint solution vector from adjoint%v_in to adjoint%v
            call chidg%data%sdata%adjoint%process_adjoint_solution(istep)

        end do


        ! Reset fields so that Rx and Jx can be computed correctly.
        ! That is, remove adjoint fields.
        call chidg%data%reset_fields()
        
        ! Compute mesh sensitivities
        call chidg%reporter('before')
        call chidg%run_adjointx(write_final=.true.)
        call chidg%reporter('after')

        ! Close ChiDG
        call chidg%shut_down('core')
        call chidg%shut_down('mpi')

        ! Delete tmp files
        if (IRANK == GLOBAL_MASTER) call delete_file(tmp_flow_file)
        if (IRANK == GLOBAL_MASTER) call delete_file(tmp_adj_file)

    end subroutine chidg_adjointx
    !******************************************************************************************





end module mod_chidg_adjointx
