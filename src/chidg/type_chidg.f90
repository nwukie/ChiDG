module type_chidg
#include <messenger.h>
    use mod_constants,              only: NFACES
    use mod_equations,              only: register_equation_builders
    use mod_operators,              only: register_operators
    use mod_models,                 only: register_models
    use mod_bc,                     only: register_bcs
    use mod_function,               only: register_functions
    use mod_grid,                   only: initialize_grid
    use mod_string,                 only: get_file_extension, string_t, get_file_prefix

    use type_chidg_data,            only: chidg_data_t
    use type_time_integrator,       only: time_integrator_t
    use mod_time,                   only: time_manager_global
    use type_linear_solver,         only: linear_solver_t
    use type_nonlinear_solver,      only: nonlinear_solver_t
    use type_preconditioner,        only: preconditioner_t
    use type_meshdata,              only: meshdata_t
    use type_bc_patch_data,         only: bc_patch_data_t
    use type_bc_state_group,        only: bc_state_group_t
    use type_bc_state,              only: bc_state_t
    use type_dict,                  only: dict_t
    use type_domain_connectivity,   only: domain_connectivity_t
    use type_partition,             only: partition_t

    use mod_time_integrators,       only: create_time_integrator
    use mod_linear_solver,          only: create_linear_solver
    use mod_nonlinear_solver,       only: create_nonlinear_solver
    use mod_preconditioner,         only: create_preconditioner

    use mod_communication,          only: establish_neighbor_communication, &
                                          establish_chimera_communication
    use mod_chidg_mpi,              only: chidg_mpi_init, chidg_mpi_finalize,   &
                                          IRANK, NRANK, ChiDG_COMM

    use mod_hdfio,                  only: read_domains_hdf, read_boundaryconditions_hdf,   &
                                          read_solution_hdf, write_solution_hdf,        &
                                          read_connectivity_hdf, read_weights_hdf,      &
                                          write_domains_hdf, read_equations_hdf
    use mod_hdf_utilities,          only: close_hdf
    use mod_partitioners,           only: partition_connectivity, send_partitions, &
                                          recv_partition
    use mpi_f08
    use mod_io
    implicit none









    !>  The ChiDG Environment container
    !!
    !!  Contains: 
    !!      - data: mesh, bcs, equations for each domain.
    !!      - a time integrator
    !!      - a nonlinear solver
    !!      - a linear solver
    !!      - a preconditioner
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------------------------
    type, public :: chidg_t

        ! Auxiliary ChiDG environment that can be used to solve sub-problems
        type(chidg_t), pointer :: auxiliary_environment

        ! Number of terms in 3D/1D solution basis expansion
        integer(ik)     :: nterms_s     = 0
        integer(ik)     :: nterms_s_1d  = 0

        ! ChiDG Files
        !type(chidg_file_t)        :: grid_file
        !type(chidg_file_t)        :: solution_file_in
        !type(chidg_file_t)        :: solution_file_out

        ! Primary data container. Mesh, equations, bc's, vectors/matrices
        type(chidg_data_t)                          :: data
        
        ! Primary algorithms, selected at run-time
        class(time_integrator_t),   allocatable     :: time_integrator
        class(nonlinear_solver_t),  allocatable     :: nonlinear_solver
        class(linear_solver_t),     allocatable     :: linear_solver
        class(preconditioner_t),    allocatable     :: preconditioner


        ! Partition of the global problem that is owned by the present ChiDG instance.
        type(partition_t)                           :: partition

        logical :: envInitialized = .false.

    contains

        ! Open/Close
        procedure   :: start_up
        procedure   :: shut_down

        ! Initialization
        procedure   :: set
        procedure   :: init

        ! Run
        procedure   :: run
        procedure   :: report

        ! IO
        procedure            :: read_grid
        procedure            :: read_domains
        procedure            :: read_boundary_conditions
        procedure            :: read_solution
        procedure            :: write_grid
        procedure            :: write_solution


    end type chidg_t
    !*******************************************************************************************





contains





    !>  ChiDG Start-Up Activities.
    !!
    !!  activity:
    !!      - 'mpi'         :: Call MPI initialization for ChiDG. This would not be called in 
    !!                         a test, since pFUnit is calling init.
    !!      - 'core'        :: Start-up ChiDG framework. Register functions, equations, 
    !!                         operators, bcs, etc.
    !!      - 'namelist'    :: Start-up Namelist IO
    !!
    !!  @author Nathan A. Wukie
    !!  @date   11/18/2016
    !!
    !!
    !-------------------------------------------------------------------------------------------
    subroutine start_up(self,activity,comm)
        class(chidg_t), intent(inout)           :: self
        character(*),   intent(in)              :: activity
        type(mpi_comm), intent(in), optional    :: comm

        integer(ik) :: ierr

        select case (trim(activity))

            !
            ! Start up MPI
            !
            case ('mpi')
                call chidg_mpi_init()


            !
            ! Start up ChiDG core
            !
            case ('core')

                ! Default communicator for 'communication' is MPI_COMM_WORLD
                if ( present(comm) ) then
                    ChiDG_COMM = comm
                else
                    ChiDG_COMM = MPI_COMM_WORLD
                end if

                ! Call environment initialization routines by default on first init call
                if (.not. self%envInitialized ) then
                    call log_init()

                ! Call environment initialization routines by default on first init call
                    ! Order matters here. Functions need to come first. Used by 
                    ! equations and bcs.
                    call register_functions()
                    call register_models()
                    call register_equation_builders()
                    call register_operators()
                    call register_bcs()
                    call initialize_grid()
                    self%envInitialized = .true.

                end if

                ! Allocate an auxiliary ChiDG environment if not already done
                if (.not. associated(self%auxiliary_environment)) then
                    allocate(self%auxiliary_environment, stat=ierr)
                    if (ierr /= 0) call AllocationError
                end if

                call self%data%time_manager%init()

                !
                ! Initialize global time_manager variable
                !
                call time_manager_global%init()


            !
            ! Start up Namelist
            !
            case ('namelist')
                call read_input()

            case default
                call chidg_signal_one(WARN,'chidg%start_up: Invalid start-up string.',trim(activity))

        end select





    end subroutine start_up
    !******************************************************************************************







    !>  Any activities that need performed before the program completely terminates.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine shut_down(self,selection)
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

                    if (allocated(self%time_integrator))  deallocate(self%time_integrator)
                    if (allocated(self%preconditioner))   deallocate(self%preconditioner)
                    if (allocated(self%linear_solver))    deallocate(self%linear_solver)
                    if (allocated(self%nonlinear_solver)) deallocate(self%nonlinear_solver)
                    call self%data%release()


                case default
                    call chidg_signal(FATAL,"chidg%shut_down: invalid shut_down string")
            end select


        else

            call log_finalize()
            call chidg_mpi_finalize()

        end if


    end subroutine shut_down
    !*****************************************************************************************






    !>  ChiDG initialization activities
    !!
    !!  activity:
    !!      - 'communication'   :: Establish local and parallel communication
    !!      - 'chimera'         :: Establish Chimera communication
    !!      - 'finalize'        :: 
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  activity   Initialization activity specification.
    !!
    !-----------------------------------------------------------------------------------------
    recursive subroutine init(self,activity)
        class(chidg_t), intent(inout)           :: self
        character(*),   intent(in)              :: activity

        character(:),   allocatable :: user_msg

        select case (trim(activity))

            !
            ! Call all initialization routines.
            !
            case ('all')
                ! geometry
                call self%init('domains')
                call self%init('bc')

                ! communication
                call self%init('comm - interior')
                call self%init('comm - chimera')

                ! matrix/vector
                call self%init('storage')
                !call self%init('finalize')



            !
            ! Initialize domain data that depend on the solution expansion
            !
            case ('domains')

                user_msg = "chidg%init('domains'): It appears the 'Solution Order' was &
                            not set for the current ChiDG instance. Try calling &
                            'call chidg%set('Solution Order',integer_input=my_order)' &
                            where my_order=1-7 indicates the solution order-of-accuracy."
                if (self%nterms_s == 0) call chidg_signal(FATAL,user_msg)

                call self%data%initialize_solution_domains(self%nterms_s)

            case ('bc')
                call self%data%initialize_solution_bc()


            !
            ! Initialize communication. Local face communication. Global parallel communication.
            !
            case ('comm - interior')
                call establish_neighbor_communication(self%data%mesh,ChiDG_COMM)


            !
            ! Initialize chimera
            !
            case ('comm - chimera')
                call establish_chimera_communication(self%data%mesh,ChiDG_COMM)


            !
            ! Initialize solver storage initialization: vectors, matrices, etc.
            !
            case ('storage')
                call self%data%initialize_solution_solver()


            !
            ! Allocate components, based on input or default input data
            !
            case ('algorithms')

                !
                ! Test chidg necessary components have been allocated
                !
                if (.not. allocated(self%time_integrator))  call chidg_signal(FATAL,"chidg%time_integrator component was not allocated")
                if (.not. allocated(self%nonlinear_solver)) call chidg_signal(FATAL,"chidg%nonlinear_solver component was not allocated")
                if (.not. allocated(self%linear_solver))    call chidg_signal(FATAL,"chidg%linear_solver component was not allocated")
                if (.not. allocated(self%preconditioner))   call chidg_signal(FATAL,"chidg%preconditioner component was not allocated")

                !
                ! Initialize preconditioner
                !
                call write_line("Initialize: preconditioner...", io_proc=GLOBAL_MASTER)
                call self%preconditioner%init(self%data)
                
                !
                ! Initialize time_integrator
                !
                call write_line("Initialize: time integrator...", io_proc=GLOBAL_MASTER)
                call self%time_integrator%init(self%data)


            case default
                call chidg_signal_one(FATAL,'chidg%init: Invalid initialization string',trim(activity))

        end select


    end subroutine init
    !*****************************************************************************************








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
    !!  @param[in]    selector    Character string for selecting the chidg component for 
    !!                            initialization
    !!  @param[in]    selection   Character string for specializing the component being 
    !!                            initialized
    !!  @param[inout] options     Dictionary for initialization options
    !!
    !-----------------------------------------------------------------------------------------
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

            case ('solution order', 'Solution Order', 'solution_order', 'Solution_Order')

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
                    !call create_time_integrator(algorithm,self%time_integrator,options)
                    call create_time_integrator(algorithm,self%time_integrator)
                else
                    !call create_time_integrator(algorithm,self%time_integrator,options)
                    call create_time_integrator(algorithm,self%time_integrator)
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
        
                

            case default
                user_msg = "chidg%set: component string was not recognized. Check spelling and that the component &
                            was registered as an option in the chidg%set routine"
                call chidg_signal_one(FATAL,user_msg,selector)


        end select


    end subroutine set
    !******************************************************************************************








    !>  Read grid from file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  gridfile        String containing a grid file name, including extension.
    !!
    !!  @param[in]  equation_set    Optionally, override the equation set for all domains
    !!  @param[in]  spacedim        Optionally, set number of spatial dimensions
    !!  
    !!  @param[in]  bc_wall         Optionally, override wall boundary functions
    !!  @param[in]  bc_inlet        Optionally, override inlet boundary functions
    !!  @param[in]  bc_outlet       Optionally, override outlet boundary functions
    !!  @param[in]  bc_symmetry     Optionally, override symmetry boundary functions
    !!  @param[in]  bc_farfield     Optionally, override farfield boundary functions
    !!  @param[in]  bc_periodic     Optionally, override periodic boundary functions
    !!
    !!  An example where overriding boundary condition is useful is computing wall distance
    !!  for a RANS calculation. First a PDE is solved using a poisson-like equation for 
    !!  wall distance and the boundary functions on walls get overridden from RANS
    !!  boundary functions to be scalar dirichlet boundary conditions. All other 
    !!  boundar functions are overridden with neumann boundary conditions.
    !!
    !------------------------------------------------------------------------------------------
    subroutine read_grid(self,gridfile,spacedim,equation_set, bc_wall, bc_inlet, bc_outlet, bc_symmetry, bc_farfield, bc_periodic, partitions_in)
        class(chidg_t),     intent(inout)               :: self
        character(*),       intent(in)                  :: gridfile
        character(*),       intent(in),     optional    :: equation_set
        integer(ik),        intent(in),     optional    :: spacedim
        class(bc_state_t),  intent(in),     optional    :: bc_wall
        class(bc_state_t),  intent(in),     optional    :: bc_inlet
        class(bc_state_t),  intent(in),     optional    :: bc_outlet
        class(bc_state_t),  intent(in),     optional    :: bc_symmetry
        class(bc_state_t),  intent(in),     optional    :: bc_farfield
        class(bc_state_t),  intent(in),     optional    :: bc_periodic
        type(partition_t),  intent(in),     optional    :: partitions_in(:)


        call write_line(' ', ltrim=.false., io_proc=GLOBAL_MASTER)
        call write_line('Reading grid... ', io_proc=GLOBAL_MASTER)


        !
        ! Read domain geometry. Also performs partitioning.
        !
        call self%read_domains(gridfile,spacedim,equation_set,partitions_in)





        !
        ! Read boundary conditions.
        !
        call self%read_boundary_conditions(gridfile, bc_wall,        &
                                                     bc_inlet,       &
                                                     bc_outlet,      &
                                                     bc_symmetry,    &
                                                     bc_farfield,    &
                                                     bc_periodic )



        !
        ! Initialize data
        !
        call self%init('all')


        call write_line('Done reading grid.', io_proc=GLOBAL_MASTER)
        call write_line(' ', ltrim=.false.,   io_proc=GLOBAL_MASTER)


    end subroutine read_grid
    !*******************************************************************************************








    !>  Read grid from file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  gridfile        String containing a grid file name, including extension.
    !!  @param[in]  spacedim        Number of spatial dimensions
    !!  @param[in]  equation_set    Optionally, specify the equation set to be initialized 
    !!                              instead of
    !!
    !!  TODO: Generalize spacedim
    !!
    !------------------------------------------------------------------------------------------
    subroutine read_domains(self,gridfile,spacedim,equation_set, partitions_in)
        class(chidg_t),     intent(inout)               :: self
        character(*),       intent(in)                  :: gridfile
        integer(ik),        intent(in),     optional    :: spacedim
        character(*),       intent(in),     optional    :: equation_set
        type(partition_t),  intent(in),     optional    :: partitions_in(:)


        type(domain_connectivity_t),    allocatable                 :: connectivities(:)
        real(rk),                       allocatable                 :: weights(:)
        type(partition_t),              allocatable, asynchronous   :: partitions(:)

        character(:),       allocatable     :: domain_equation_set
        type(meshdata_t),   allocatable     :: meshdata(:)
        integer(ik)                         :: idom, iread, ierr, &
                                               domain_dimensionality, ielem, eqn_ID


        call write_line(' ',                      ltrim=.false., io_proc=GLOBAL_MASTER)
        call write_line('   Reading domains... ', ltrim=.false., io_proc=GLOBAL_MASTER)



        !
        ! Partitions defined from user input
        !
        if ( present(partitions_in) ) then

            self%partition = partitions_in(IRANK+1)


        !
        ! Partitions from partitioner tool
        !
        else

            !
            ! Master rank: Read connectivity, partition connectivity, distribute partitions
            !
            call write_line("   partitioning...", ltrim=.false., io_proc=GLOBAL_MASTER)
            if ( IRANK == GLOBAL_MASTER ) then

                call read_connectivity_hdf(gridfile,connectivities)
                call read_weights_hdf(gridfile,weights)

                call partition_connectivity(connectivities, weights, partitions)

                call send_partitions(partitions,MPI_COMM_WORLD)
            end if


            !
            ! All ranks: Receive partition from GLOBAL_MASTER
            !
            call write_line("   distributing partitions...", ltrim=.false., io_proc=GLOBAL_MASTER)
            call recv_partition(self%partition,MPI_COMM_WORLD)


        end if ! partitions in from user



        !
        ! Read data from hdf file
        !
        do iread = 0,NRANK-1
            if ( iread == IRANK ) then

                call read_equations_hdf(gridfile, self%data)
                call read_domains_hdf(gridfile,self%partition,meshdata)

            end if
            call MPI_Barrier(ChiDG_COMM,ierr)
        end do



        !
        ! Add domains to ChiDG%data
        !
        call write_line("   processing...", ltrim=.false., io_proc=GLOBAL_MASTER)
        do idom = 1,size(meshdata)


            ! Use spacedim if specified, else default to 3D
            if (present(spacedim)) then
                domain_dimensionality = spacedim
            else 
                domain_dimensionality = 3
            end if


            ! Use equation_set if specified, else default to the grid file data
            if (present(equation_set)) then
                domain_equation_set = equation_set
                call self%data%add_equation_set(equation_set)
            else
                domain_equation_set = meshdata(idom)%eqnset
            end if



            ! Get the equation set identifier
            eqn_ID = self%data%get_equation_set_id(domain_equation_set)

            call self%data%mesh%add_domain( trim(meshdata(idom)%name),    &
                                            meshdata(idom)%points,        &
                                            meshdata(idom)%connectivity,  &
                                            meshdata(idom)%nelements_g,   &
                                            domain_dimensionality,        &
                                            meshdata(idom)%nterms_c,      &
                                            meshdata(idom)%coord_system,  &
                                            eqn_ID )



        end do !idom



        call write_line('   Done reading domains... ', ltrim=.false., io_proc=GLOBAL_MASTER)
        call write_line(' ',                           ltrim=.false., io_proc=GLOBAL_MASTER)

    end subroutine read_domains
    !******************************************************************************************











    !>  Read boundary conditions from grid file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!  @param[in]  gridfile    String specifying a gridfile, including extension.
    !!
    !-----------------------------------------------------------------------------------------
    subroutine read_boundary_conditions(self, gridfile, bc_wall, bc_inlet, bc_outlet, bc_symmetry, bc_farfield, bc_periodic)
        class(chidg_t),     intent(inout)               :: self
        character(*),       intent(in)                  :: gridfile
        class(bc_state_t),  intent(in),     optional    :: bc_wall
        class(bc_state_t),  intent(in),     optional    :: bc_inlet
        class(bc_state_t),  intent(in),     optional    :: bc_outlet
        class(bc_state_t),  intent(in),     optional    :: bc_symmetry
        class(bc_state_t),  intent(in),     optional    :: bc_farfield
        class(bc_state_t),  intent(in),     optional    :: bc_periodic

        type(bc_patch_data_t),  allocatable     :: bc_patch_data(:)
        type(string_t)                          :: bc_group_name
        type(bc_state_group_t), allocatable     :: bc_state_groups(:)
        type(string_t)                          :: group_name
        integer(ik)                             :: idom, ndomains, iface, ibc, ierr, iread, bc_ID


        call write_line(' ',                                  ltrim=.false., io_proc=GLOBAL_MASTER)
        call write_line('   Reading boundary conditions... ', ltrim=.false., io_proc=GLOBAL_MASTER)


        !
        ! Call boundary condition reader based on file extension
        !
        do iread = 0,NRANK-1
            if ( iread == IRANK ) then

                call read_boundaryconditions_hdf(gridfile,bc_patch_data,bc_state_groups,self%partition)

            end if
            call MPI_Barrier(ChiDG_COMM,ierr)
        end do





        !
        ! Add all boundary condition state groups
        !
        call write_line('   processing groups...', ltrim=.false., io_proc=GLOBAL_MASTER)
        do ibc = 1,size(bc_state_groups)

            call self%data%add_bc_state_group(bc_state_groups(ibc), bc_wall     = bc_wall,      &
                                                                    bc_inlet    = bc_inlet,     &
                                                                    bc_outlet   = bc_outlet,    &
                                                                    bc_symmetry = bc_symmetry,  &
                                                                    bc_farfield = bc_farfield,  &
                                                                    bc_periodic = bc_periodic )
      
        end do !ibc



        !
        ! Add boundary condition patch groups
        !
        call write_line('   processing patches...', ltrim=.false., io_proc=GLOBAL_MASTER)
        ndomains = size(bc_patch_data)
        do idom = 1,ndomains
            do iface = 1,NFACES

                bc_group_name = bc_patch_data(idom)%bc_group_name%at(iface)
                bc_ID         = self%data%get_bc_state_group_id(bc_group_name%get())

                call self%data%mesh%add_bc_patch(bc_patch_data(idom)%domain_name,               &
                                                 bc_group_name%get(),                           &
                                                 bc_patch_data(idom)%bc_connectivity(iface),    &
                                                 bc_ID)

            end do !iface
        end do !ipatch



        call write_line('   Done reading boundary conditions.', ltrim=.false., io_proc=GLOBAL_MASTER)
        call write_line(' ',                                    ltrim=.false., io_proc=GLOBAL_MASTER)


    end subroutine read_boundary_conditions
    !*****************************************************************************************











    !>  Read solution from file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  solutionfile    String containing a solution file name, including extension.
    !!
    !-----------------------------------------------------------------------------------------
    subroutine read_solution(self,file_name)
        class(chidg_t),     intent(inout)           :: self
        character(*),       intent(in)              :: file_name

        integer(ik) :: iread, ierr

        call write_line(' ', ltrim=.false.,      io_proc=GLOBAL_MASTER)
        call write_line(' Reading solution... ', io_proc=GLOBAL_MASTER)



        !
        ! Read solution from hdf file
        !
        call write_line("   reading from: ", file_name, ltrim=.false., io_proc=GLOBAL_MASTER)

        call read_solution_hdf(file_name,self%data)


        call write_line('Done reading solution.', io_proc=GLOBAL_MASTER)
        call write_line(' ', ltrim=.false.,       io_proc=GLOBAL_MASTER)

    end subroutine read_solution
    !*****************************************************************************************










!    !>  Check the equation_set's for any auxiliary fields that are required. 
!    !!
!    !!      #1: Check equation_set's for auxiliary fields
!    !!      #2: For auxiliary fields that are found, check the file for the auxiliary field.
!    !!      #3: If no auxiliary field in file, check auxiliary drivers for a rule to compute the field
!    !!      #4: If no rule, error. We don't have the field, and we also don't know how to compute it.
!    !!
!    !!  @author Nathan A. Wukie
!    !!  @date   11/21/2016
!    !!
!    !!
!    !!
!    !-----------------------------------------------------------------------------------------
!    subroutine check_auxiliary_fields(self)
!        class(chidg_t), intent(inout)   :: self
!
!        type(string_t), allocatable :: auxiliary_fields(:)
!        logical,        allocatable :: domain_needs_aux_field(:)
!        logical,        allocatable :: file_has_aux_field(:)
!        logical,        allocatable :: field_in_file
!
!
!
!!        aux_fields = self%data%get_auxiliary_fields()
!!
!!
!!        do iaux = 1,size(aux_fields)
!!
!!
!!            !
!!            ! Check which domains use the auxiliary field
!!            !
!!            do idom = 1,self%data%ndomains()
!!                domain_uses_field(idom) = self%data%eqnset(idom)%uses_auxiliary_field(aux_fields(iaux))
!!            end do !idom
!!
!!
!!            !
!!            do idom = 1,self%data%ndomains()
!!                file_has_aux_field(idom) = file%domain_has_field(idom,aux_field(iaux))
!!            end do !idom
!!
!!
!!            !
!!            ! If any domain doesn't have the field in file, get from a pre-defined rule
!!            !
!!            if (any(file_has_aux_field == .false.)) then
!!                call initialize_auxiliary_field(aux_field(iaux))
!!            end if
!!
!!
!!        end do !iaux
!
!        
!
!
!    end subroutine check_auxiliary_fields
!    !*****************************************************************************************







    !>  Write grid to file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   3/21/2016
    !!
    !!  @param[in]  file    String containing a solution file name, including extension.
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine write_grid(self,file_name)
        class(chidg_t),     intent(inout)           :: self
        character(*),       intent(in)              :: file_name


        call write_line(' ', ltrim=.false., io_proc=GLOBAL_MASTER)
        call write_line('Writing grid... ', io_proc=GLOBAL_MASTER)


        !
        ! Call grid reader based on file extension
        !
        call write_line("   writing to: ", file_name, ltrim=.false., io_proc=GLOBAL_MASTER)
        call write_domains_hdf(self%data,file_name)


        call write_line("Done writing grid.", io_proc=GLOBAL_MASTER)
        call write_line(' ', ltrim=.false.,   io_proc=GLOBAL_MASTER)

    end subroutine write_grid
    !*****************************************************************************************









    !>  Write solution to file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]  solutionfile    String containing a solution file name, including extension.
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine write_solution(self,file_name)
        class(chidg_t),     intent(inout)           :: self
        character(*),       intent(in)              :: file_name


        call write_line(' ', ltrim=.false.,     io_proc=GLOBAL_MASTER)
        call write_line('Writing solution... ', io_proc=GLOBAL_MASTER)


        !
        ! Call grid reader based on file extension
        !
        call write_line("   writing to:", file_name, ltrim=.false., io_proc=GLOBAL_MASTER)

        call write_solution_hdf(self%data,file_name)
        call self%time_integrator%write_time_options(self%data,file_name)



        call write_line("Done writing solution.", io_proc=GLOBAL_MASTER)
        call write_line(' ', ltrim=.false.,       io_proc=GLOBAL_MASTER)

    end subroutine write_solution
    !*****************************************************************************************












    !>  Run ChiDG simulation
    !!
    !!      - This routine passes the domain data, nonlinear solver, linear solver, and 
    !!        preconditioner components to the time integrator to take a step.
    !!
    !!  Optional input parameters:
    !!      write_initial   control writing initial solution to file. Default: .false.
    !!      write_final     control writing final solution to file. Default: .true.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!  @date   2/7/2017
    !!
    !------------------------------------------------------------------------------------------
    subroutine run(self, write_initial, write_final)
        class(chidg_t), intent(inout)           :: self
        logical,        intent(in), optional    :: write_initial
        logical,        intent(in), optional    :: write_final

        character(100)              :: filename
        character(:),   allocatable :: prefix
        integer(ik)                 :: istep, nsteps, wcount
        logical                     :: option_write_initial, option_write_final

    

        call write_line("---------------------------------------------------", io_proc=GLOBAL_MASTER)
        call write_line("                                                   ", io_proc=GLOBAL_MASTER, delimiter='none')
        call write_line("           Running ChiDG simulation...             ", io_proc=GLOBAL_MASTER, delimiter='none')
        call write_line("                                                   ", io_proc=GLOBAL_MASTER, delimiter='none')
        call write_line("---------------------------------------------------", io_proc=GLOBAL_MASTER)


        !
        ! Initialize algorithms
        !
        call self%init('algorithms')



        !
        ! Check optional incoming parameters
        !
        if (present(write_initial)) then
            option_write_initial = write_initial
        else
            option_write_initial = .false.
        end if

        if (present(write_final)) then
            option_write_final = write_final
        else
            option_write_final = .true.
        end if




        !
        ! Write initial solution
        !
        if (option_write_initial) call self%write_solution('initial.h5')




        
        !
        ! Initialize time integrator state
        !
        call self%time_integrator%initialize_state(self%data)

        
        !
        ! Get the prefix in the file name in case of multiple output files
        !
        prefix = get_file_prefix(solutionfile_out,'.h5')       




        !
        ! Execute time_integrator, nsteps times 
        !
        wcount = 1
        nsteps = self%data%time_manager%nsteps
        do istep = 1,nsteps
            

            call write_line("- Step ", istep, io_proc=GLOBAL_MASTER)


            !
            ! 1: Update time t
            ! 2: Call time integrator to take a step
            !
            self%data%time_manager%t = self%data%time_manager%dt*istep
            call self%time_integrator%step(self%data,self%nonlinear_solver, &
                                                     self%linear_solver,    &
                                                     self%preconditioner)
           
           

            !
            ! Write solution every nwrite steps
            !
            if (wcount == self%data%time_manager%nwrite) then
                write(filename, "(A,I7.7,A3)") trim(prefix)//'_', istep, '.h5'
                call self%write_solution(filename)
                wcount = 0
            end if



            !
            ! Print diagnostics
            !
            call write_line("-  System residual |R(Q)|: ", self%time_integrator%residual_norm%at(istep), delimiter='', io_proc=GLOBAL_MASTER)


            wcount = wcount + 1


        end do !istep

                
        !
        ! Write the final solution to hdf file
        !        
        if (option_write_final) call self%write_solution(solutionfile_out)


    end subroutine run
    !*****************************************************************************************










    !> Report on ChiDG simulation
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!
    !!
    !------------------------------------------------------------------------------------------
    subroutine report(self,selection)
        class(chidg_t), intent(inout)   :: self
        character(*),   intent(in)      :: selection

        integer(ik) :: ireport, ierr


        if ( trim(selection) == 'before' ) then


            do ireport = 0,NRANK-1
                if ( ireport == IRANK ) then
                    call write_line('MPI Rank: ', IRANK, io_proc=IRANK)
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
    !*****************************************************************************************



















end module type_chidg
