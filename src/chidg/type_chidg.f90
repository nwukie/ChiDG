module type_chidg
#include <messenger.h>
    use mod_constants,              only: NFACES
    use mod_equations,              only: register_equation_builders
    use mod_operators,              only: register_operators
    use mod_bc,                     only: register_bcs
    use mod_function,               only: register_functions
    use mod_grid,                   only: initialize_grid
    use mod_io,                     only: read_input
    use mod_string,                 only: get_file_extension, string_t

    use type_chidg_data,            only: chidg_data_t
    use type_time_integrator,       only: time_integrator_t
    use type_linear_solver,         only: linear_solver_t
    use type_nonlinear_solver,      only: nonlinear_solver_t
    use type_preconditioner,        only: preconditioner_t
    use type_meshdata,              only: meshdata_t
    use type_bc_patch_data,         only: bc_patch_data_t
    use type_bc_group_data,         only: bc_group_data_t
    use type_bc_state,              only: bc_state_t
    use type_dict,                  only: dict_t
    use type_domain_connectivity,   only: domain_connectivity_t
    use type_partition,             only: partition_t

    use mod_time_integrators,       only: create_time_integrator
    use mod_linear_solver,          only: create_linear_solver
    use mod_nonlinear_solver,       only: create_nonlinear_solver
    use mod_preconditioner,         only: create_preconditioner

    use mod_communication,          only: establish_neighbor_communication, establish_chimera_communication
    use mod_chidg_mpi,              only: chidg_mpi_init, chidg_mpi_finalize, ChiDG_COMM, IRANK, NRANK
    use mpi_f08

    use mod_hdfio,                  only: read_grid_hdf, read_boundaryconditions_hdf, &
                                          read_solution_hdf, write_solution_hdf, read_connectivity_hdf, &
                                          read_weights_hdf
    use mod_hdf_utilities,          only: close_hdf
    use mod_partitioners,           only: partition_connectivity, send_partitions, recv_partition
    implicit none









    !>  The ChiDG Environment container
    !!
    !!      - Contains an array of domains, a time integrator, a nonlinear solver, a linear solver, and a preconditioner
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !--------------------------------------------------------------------------------------------------------
    type, public    :: chidg_t


        integer(ik) :: ntime        = 1 ! Number of time instances being solved for
        integer(ik) :: nterms_s     = 0 ! Number of terms in the 3D solution basis expansion
        integer(ik) :: nterms_s_1d  = 0 ! Number of terms in the 1D solution basis expansion


        ! Primary data container. Mesh, equations, bc's, vectors/matrices
        type(chidg_data_t)   :: data
        
        ! Primary algorithms, selected at run-time
        class(time_integrator_t),   allocatable     :: time_integrator
        class(nonlinear_solver_t),  allocatable     :: nonlinear_solver
        class(linear_solver_t),     allocatable     :: linear_solver
        class(preconditioner_t),    allocatable     :: preconditioner


        ! Partition of the global problem that is owned by the present ChiDG instance.
        type(partition_t)                           :: partition

        logical :: envInitialized = .false.

    contains


        procedure   :: init
        procedure   :: close
        procedure   :: set

        procedure   :: run
        procedure   :: report

        ! IO procedures
        procedure   :: read_grid
        procedure   :: read_boundaryconditions
        procedure   :: read_solution
        procedure   :: write_solution

        ! Initialization
        procedure   :: initialize_solution_domains
        procedure   :: initialize_solution_solver

    end type chidg_t
    !*********************************************************************************************************







contains



    !> ChiDG environment initialization routine
    !!      - Call initiailization procedures for equations, grid data, reading input
    !!      - chidg%init('env') should be called before any activity with ChiDG is begun.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  level   Initialization level specification. 'env', 'communication', 'io', or 'finalize'
    !!
    !--------------------------------------------------------------------------------------------
    subroutine init(self,level,comm)
        class(chidg_t), intent(inout)           :: self
        character(*),   intent(in)              :: level
        type(mpi_comm), intent(in), optional    :: comm




        ! Valid strings are:
        !   - 'env'             Basic environment initialization. Equations and supporting grid data
        !   - 'mpi'             Call MPI initialization for ChiDG. This would not be called in a test, since pFUnit is calling init.
        !   - 'communication'   Establish local and global communication
        !   - 'io'              Read namelist input.
        !   - 'finalize'        Call component initialization routines before run

        ! Call environment initialization routines by default on first init call
        if (.not. self%envInitialized ) then
            call log_init()


            ! Order matters here. Functions need to come first. Used by equations and bcs.
            call register_functions()
            call register_equation_builders()
            call register_operators()
            call register_bcs()

            call initialize_grid()

            self%envInitialized = .true.

        end if





        select case (trim(level))

            !
            ! Do nothing else. Just to ensure the environment commands above are executed
            !
            case ('env')

                !
                ! Default communicator for 'communication' is MPI_COMM_WORLD
                !
                if ( present(comm) ) then
                    ChiDG_COMM = comm
                else
                    ChiDG_COMM = MPI_COMM_WORLD
                end if


            case ('mpi')
                call chidg_mpi_init()

            !
            ! Read Namelist input file
            !
            case ('io')
                call read_input()


            !
            ! Initialize communication. Local face communication. Global parallel communication.
            !
            case ('communication')
                call establish_neighbor_communication(self%data%mesh,ChiDG_COMM)


            !
            ! Initialize chimera
            !
            case ('chimera')
                call establish_chimera_communication(self%data%mesh,ChiDG_COMM)


            !
            ! Allocate components, based on input or default input data
            !
            case ('finalize')

                !
                ! Test chidg necessary components have been allocated
                !
                if (.not. allocated(self%time_integrator))  call chidg_signal(FATAL,"chidg%time_integrator component was not allocated")
                if (.not. allocated(self%nonlinear_solver)) call chidg_signal(FATAL,"chidg%nonlinear_solver component was not allocated")
                if (.not. allocated(self%linear_solver))    call chidg_signal(FATAL,"chidg%linearsolver component was not allocated")
                if (.not. allocated(self%preconditioner))   call chidg_signal(FATAL,"chidg%preconditioner component was not allocated")


                call self%time_integrator%init(self%data)
                call self%preconditioner%init(self%data)



            case default
                call chidg_signal(WARN,'chidg_t: Invalid initialization string')

        end select


    end subroutine init
    !**********************************************************************************************************








    !>  Set ChiDG environment components
    !!
    !!      -   Set time-integrator
    !!      -   Set nonlinear solver
    !!      -   Set linear solver
    !!      -   Set preconditioner
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]      selector    Character string for selecting the chidg component for initialization
    !!  @param[in]      selection   Character string for specializing the component being initialized
    !!  @param[inout]   options     Dictionary for initialization options
    !!
    !--------------------------------------------------------------------------------------------------------
    subroutine set(self,selector,algorithm,integer_input,real_input,options)
        class(chidg_t),         intent(inout)   :: self
        character(*),           intent(in)      :: selector
        character(*), optional, intent(in)      :: algorithm
        integer(ik),  optional, intent(in)      :: integer_input
        real(rk),     optional, intent(in)      :: real_input
        type(dict_t), optional, intent(inout)   :: options 

        character(:),   allocatable :: user_msg
        integer(ik)                 :: ierr


        !
        ! Check options for the algorithm family of inputs
        !
        select case (trim(selector))
            
            case ('Time','time','time_integrator','Time_Integrator','timeintegrator','TimeIntegrator', 'time integrator', 'Time Integrator', &
                  'nonlinearsolver','NonlinearSolver','nonlinear_solver','Nonlinear_Solver','nonlinear solver', 'Nonlinear Solver', &
                  'linearsolver','LinearSolver','linear_solver','Linear_Solver','linear solver', 'Linear Solver', &
                  'preconditioner','Preconditioner')

                user_msg = "chidg%set: The component being set needs an algorithm string passed in &
                            along with it. Try 'call chidg%set('your component', algorithm='your algorithm string')"
                if (.not. present(algorithm)) call chidg_signal_one(FATAL,user_msg,trim(selector))

        end select


        !
        ! Check options for integer family of inputs
        !
        select case (trim(selector))

            case ('solution order', 'Solution Order', 'solution_order', 'Solution_Order', &
                  'ntime', 'Number of Time Instances', 'NTIME', 'n time')

                user_msg = "chidg%set: The component being set needs an integer passed in &
                            along with it. Try 'call chidg%set('your component', integer_input=my_int)"
                if (.not. present(integer_input)) call chidg_signal_one(FATAL,user_msg,trim(selector))

        end select






        !
        ! Actually go in and call the specialized routine based on 'selector'
        !
        select case (trim(selector))
            !
            ! Allocation for time integrator
            !
            case ('Time','time','time_integrator','Time_Integrator','timeintegrator','TimeIntegrator', 'time integrator', 'Time Integrator')
                if (allocated(self%time_integrator)) then
                    deallocate(self%time_integrator)
                    call create_time_integrator(algorithm,self%time_integrator,options)
                else
                    call create_time_integrator(algorithm,self%time_integrator,options)
                end if



            !
            ! Allocation for nonlinear solver
            !
            case ('nonlinearsolver','NonlinearSolver','nonlinear_solver','Nonlinear_Solver','nonlinear solver', 'Nonlinear Solver')
                if (allocated(self%nonlinear_solver)) then
                    deallocate(self%nonlinear_solver)
                    call create_nonlinear_solver(algorithm,self%nonlinear_solver,options)
                else
                    call create_nonlinear_solver(algorithm,self%nonlinear_solver,options)
                end if




            !
            ! Allocation for linear solver
            !
            case ('linearsolver','LinearSolver','linear_solver','Linear_Solver','linear solver', 'Linear Solver')
                if (allocated(self%linear_solver)) then
                    deallocate(self%linear_solver)
                    call create_linear_solver(algorithm,self%linear_solver,options)
                else
                    call create_linear_solver(algorithm,self%linear_solver,options)
                end if


            !
            ! Allocation for preconditioner
            !
            case ('preconditioner','Preconditioner')
                if (allocated(self%preconditioner)) deallocate(self%preconditioner)

                call create_preconditioner(algorithm,self%preconditioner)



            !
            ! Set the 'solution order'. Order-of-accuracy, that is. Compute the number of terms
            ! in the 1D, 3D solution bases
            !
            case ('solution order', 'Solution Order', 'solution_order', 'Solution_Order')
                self%nterms_s_1d = integer_input
                self%nterms_s    = self%nterms_s_1d * self%nterms_s_1d * self%nterms_s_1d
        

            case ('ntime', 'Number of Time Instances', 'NTIME', 'n time')
                self%ntime = integer_input
                

            case default
                user_msg = "chidg%set: component string was not recognized. Check spelling and that the component &
                            was registered as an option in the chidg%set routine"
                call chidg_signal_one(FATAL,user_msg,selector)


        end select


    end subroutine set
    !********************************************************************************************************








    !>  Read grid from file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  gridfile    String containing a grid file name, including extension.
    !!  @param[in]  spacedim    Number of spatial dimensions
    !!  @param[in]  equation_set    Optionally, specify the equation set to be initialized instead of
    !!
    !!  TODO: Generalize spacedim
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine read_grid(self,gridfile,spacedim,equation_set)
        class(chidg_t),             intent(inout)           :: self
        character(*),               intent(in)              :: gridfile
        character(*),   optional,   intent(in)              :: equation_set
        integer(ik),    optional,   intent(in)              :: spacedim



        type(domain_connectivity_t),    allocatable :: connectivities(:)
        real(rk),                       allocatable :: weights(:)
        type(partition_t),              allocatable :: partitions(:)

        character(len=5),   dimension(1)    :: extensions
        character(len=:),   allocatable     :: extension, domain_equation_set
        type(meshdata_t),   allocatable     :: meshdata(:)
        integer                             :: iext, extloc, idom, ndomains, iread, ierr, domain_dimensionality


        if ( IRANK == GLOBAL_MASTER ) call write_line("Reading grid")


        !
        ! Master rank: Read connectivity, partition connectivity, distribute partitions
        !
        if ( IRANK == GLOBAL_MASTER ) then

            call read_connectivity_hdf(gridfile,connectivities)
            call read_weights_hdf(gridfile,weights)

            call partition_connectivity(connectivities, weights, partitions)

            call send_partitions(partitions,MPI_COMM_WORLD)
        end if


        !
        ! All ranks: Receive partition from GLOBAL_MASTER
        !
        call recv_partition(self%partition,MPI_COMM_WORLD)







        !
        ! Get filename extension
        !
        extensions = ['.h5']
        extension = get_file_extension(gridfile, extensions)


        !
        ! Call grid reader based on file extension
        !
        do iread = 0,NRANK-1
            if ( iread == IRANK ) then


                if ( extension == '.h5' ) then 
                    call read_grid_hdf(gridfile,self%partition,meshdata)
                else
                    call chidg_signal(FATAL,"chidg%read_grid: grid file extension not recognized")
                end if


            end if
            call MPI_Barrier(ChiDG_COMM,ierr)
        end do



        !
        ! Add domains to ChiDG%data
        !
        ndomains = size(meshdata)
        do idom = 1,ndomains



            ! Use spacedim if specified, else default to 3D
            if (present(spacedim)) then
                domain_dimensionality = spacedim
            else 
                domain_dimensionality = 3
            end if


            ! Use equation_set if specified, else default to the grid file data
            if (present(equation_set)) then
                domain_equation_set = equation_set
            else
                domain_equation_set = meshdata(idom)%eqnset
            end if



            call self%data%add_domain(                              &
                                      trim(meshdata(idom)%name),    &
                                      meshdata(idom)%points,        &
                                      meshdata(idom)%connectivity,  &
                                      domain_dimensionality,        &
                                      meshdata(idom)%nterms_c,      &
                                      domain_equation_set           &
                                      )

        end do


    end subroutine read_grid
    !**********************************************************************************************************











    !>  Read boundary conditions from grid file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!  @param[in]  gridfile    String specifying a gridfile, including extension.
    !!
    !-------------------------------------------------------------------------------------------------------------
    subroutine read_boundaryconditions(self, gridfile, bc_wall, bc_inlet, bc_outlet, bc_symmetry, bc_farfield)
        class(chidg_t),     intent(inout)               :: self
        character(*),       intent(in)                  :: gridfile
        class(bc_state_t),  intent(in),     optional    :: bc_wall
        class(bc_state_t),  intent(in),     optional    :: bc_inlet
        class(bc_state_t),  intent(in),     optional    :: bc_outlet
        class(bc_state_t),  intent(in),     optional    :: bc_symmetry
        class(bc_state_t),  intent(in),     optional    :: bc_farfield

        character(len=5),   dimension(1)    :: extensions
        character(len=:),   allocatable     :: extension
        type(bc_patch_data_t),  allocatable :: bc_patches(:)
        type(bc_group_data_t),  allocatable :: bc_groups(:)
        type(string_t)                      :: group_name
        integer                             :: idom, ndomains, iface, ierr, iread

        if (IRANK == GLOBAL_MASTER) call write_line('Reading boundary conditions')

        !
        ! Get filename extension
        !
        extensions = ['.h5']
        extension = get_file_extension(gridfile, extensions)


        !
        ! Call boundary condition reader based on file extension
        !
        do iread = 0,NRANK-1
            if ( iread == IRANK ) then


                if ( extension == '.h5' ) then
                    call read_boundaryconditions_hdf(gridfile,bc_patches,bc_groups,self%partition)
                else
                    call chidg_signal(FATAL,"chidg%read_boundaryconditions: grid file extension not recognized")
                end if


            end if
            call MPI_Barrier(ChiDG_COMM,ierr)
        end do



        !
        ! Add boundary conditions to ChiDG
        !
        ndomains = size(bc_patches)
        do idom = 1,ndomains
            do iface = 1,NFACES

                group_name = bc_patches(idom)%bc_group%at(iface)
                call self%data%add_bc(bc_patches(idom)%domain_, bc_patches(idom)%bc_connectivity(iface), group_name%get(),bc_groups)
!                call self%data%add_bc(bcdata(idom)%domain_, bcdata(idom)%bcs(iface), bcdata(idom)%bc_connectivity(iface), &
!                                      bc_wall,      &
!                                      bc_inlet,     &
!                                      bc_outlet,    &
!                                      bc_symmetry,  &
!                                      bc_farfield)

            end do !iface
        end do !idom


    end subroutine read_boundaryconditions
    !**************************************************************************************************************












    !>  Read solution from file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  solutionfile    String containing a solution file name, including extension.
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine read_solution(self,solutionfile)
        class(chidg_t),     intent(inout)           :: self
        character(*),       intent(in)              :: solutionfile

        character(len=5),   dimension(1)    :: extensions
        character(len=:),   allocatable     :: extension
        type(meshdata_t),   allocatable     :: solutiondata(:)
        integer                             :: iext, extloc, idom, ndomains, iread, ierr


        !
        ! Get filename extension
        !
        extensions = ['.h5']
        extension = get_file_extension(solutionfile, extensions)


        !
        ! Call grid reader based on file extension
        !
        do iread = 0,NRANK-1
            if ( iread == IRANK ) then

                if ( extension == '.h5' ) then
                    call read_solution_hdf(solutionfile,self%data)
                else
                    call chidg_signal(FATAL,"chidg%read_solution: grid file extension not recognized")
                end if

            end if
            call MPI_Barrier(ChiDG_COMM,ierr)
        end do ! iread


    end subroutine read_solution
    !************************************************************************************************************












    

    !>  Initialize all solution and solver storage.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2016
    !!
    !!  @param[in]  nterms_s    Number of terms in the solution polynomial expansion.
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine initialize_solution_domains(self)
        class(chidg_t),     intent(inout)   :: self
!        integer(ik),        intent(in)      :: nterms_s

        character(:),   allocatable :: user_msg

        !
        ! TODO: put in checks for prerequisites
        !


        !
        ! Check that the order for the solution basis expansion has been set.
        !
        user_msg = "chidg%initialize_solution_domains: It appears the 'Solution Order' was not set for the &
                    current ChiDG instance. Try calling 'call chidg%set('Solution Order', integer_input=my_order)' &
                    where my_order=1-7 indicates the solution order-of-accuracy."
        if (self%nterms_s == 0) call chidg_signal(FATAL,user_msg)

        !
        ! Call domain solution storage initialization: mesh data structures that depend on solution expansion etc.
        !
        call self%data%initialize_solution_domains(self%nterms_s,self%ntime)


    end subroutine initialize_solution_domains
    !************************************************************************************************************






    !>  Initialize all solution and solver storage.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   4/11/2016
    !!
    !!  @param[in]  nterms_s    Number of terms in the solution polynomial expansion.
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine initialize_solution_solver(self)
        class(chidg_t),     intent(inout)   :: self


        !
        ! TODO: put in checks for prerequisites
        !

        !
        ! Call solver solution storage initialization: vectors, matrices, etc.
        !
        call self%data%initialize_solution_solver()



    end subroutine initialize_solution_solver
    !************************************************************************************************************
















    !> Write solution to file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  solutionfile    String containing a solution file name, including extension.
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine write_solution(self,solutionfile)
        class(chidg_t),     intent(inout)           :: self
        character(*),       intent(in)              :: solutionfile

        character(len=5),   dimension(1)    :: extensions
        character(len=:),   allocatable     :: extension
        type(meshdata_t),   allocatable     :: solutiondata(:)
        integer                             :: iext, extloc, idom, ndomains, iwrite, ierr


        !
        ! Get filename extension
        !
        extensions = ['.h5']
        extension = get_file_extension(solutionfile, extensions)


        !
        ! Call grid reader based on file extension
        !
        do iwrite = 0,NRANK-1
            if ( iwrite == IRANK ) then


                if ( extension == '.h5' ) then
                    call write_solution_hdf(solutionfile,self%data)
                else
                    call chidg_signal(FATAL,"chidg%write_solution: grid file extension not recognized")
                end if


            end if
            call MPI_Barrier(ChiDG_COMM,ierr)
        end do


    end subroutine write_solution
    !************************************************************************************************************














    !>  Run ChiDG simulation
    !!
    !!      - This routine passes the domain data, nonlinear solver, linear solver, and preconditioner
    !!        components to the time integrator for iteration
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine run(self)
        class(chidg_t),     intent(inout)   :: self


        call self%time_integrator%iterate(self%data,self%nonlinear_solver,self%linear_solver,self%preconditioner)


    end subroutine run
    !************************************************************************************************************










    !> Report on ChiDG simulation
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine report(self,selection)
        class(chidg_t), intent(inout)   :: self
        character(*),   intent(in)      :: selection

        integer(ik) :: ireport, ierr


        if ( trim(selection) == 'before' ) then


            do ireport = 0,NRANK-1
                if ( ireport == IRANK ) then
                    call write_line('MPI Rank: ', IRANK)
                    call self%data%report('grid')
                end if
                call MPI_Barrier(ChiDG_COMM,ierr)
            end do









        else if ( trim(selection) == 'after' ) then

            if ( IRANK == GLOBAL_MASTER ) then
                !call self%time_integrator%report()
                call self%nonlinear_solver%report()
                !call self%preconditioner%report()
            end if

        end if


    end subroutine report
    !************************************************************************************************************









    !>  Any activities that need performed before the program completely terminates.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine close(self,selection)
        class(chidg_t), intent(inout)               :: self
        character(*),   intent(in),     optional    :: selection



        if ( present(selection) ) then 
            select case (selection)
                case ('log')
                    call log_finalize()
                case ('mpi')
                    call chidg_mpi_finalize()
                case ('core')   ! All except mpi
                    call log_finalize()
                    call close_hdf()

                case default
                    call chidg_signal(FATAL,"chidg%close: invalid close string")
            end select


        else

            call log_finalize()
            call chidg_mpi_finalize()

        end if


    end subroutine close
    !************************************************************************************************************











end module type_chidg
