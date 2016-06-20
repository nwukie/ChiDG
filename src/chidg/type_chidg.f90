module type_chidg
#include <messenger.h>
    use mod_constants,          only: NFACES
    use mod_equations,          only: register_equations
    use mod_bc,                 only: register_bcs
    use mod_function,           only: register_functions
    use mod_grid,               only: initialize_grid
    use mod_io,                 only: read_input
    use mod_string_utilities,   only: get_file_extension

    use type_chidg_data,        only: chidg_data_t
    use type_time_scheme,       only: time_scheme_t
    use type_linear_solver,     only: linear_solver_t
    use type_nonlinear_solver,  only: nonlinear_solver_t
    use type_preconditioner,    only: preconditioner_t
    use type_meshdata,          only: meshdata_t
    use type_bcdata,            only: bcdata_t
    use type_dict,              only: dict_t
    use type_partition,         only: partition_t

    use mod_time_scheme,        only: create_time_scheme
    use mod_linear_solver,      only: create_linear_solver
    use mod_nonlinear_solver,   only: create_nonlinear_solver
    use mod_preconditioner,     only: create_preconditioner
    use mod_chimera,            only: detect_chimera_faces,  &
                                      detect_chimera_donors, &
                                      compute_chimera_interpolators
    use mod_communication,      only: establish_communication
    use mod_chidg_mpi,          only: chidg_mpi_init, chidg_mpi_finalize

    use mod_hdfio,              only: read_grid_hdf, read_grid_partition_hdf, read_boundaryconditions_hdf, read_boundaryconditions_partition_hdf, read_solution_hdf, write_solution_hdf
    implicit none





    !>  The ChiDG Environment container
    !!
    !!      - Contains an array of domains, a time advancement scheme, a nonlinear solver, a linear solver, and a preconditioner
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !--------------------------------------------------------------------------------------------------------
    type, public    :: chidg_t

        type(chidg_data_t)                          :: data
        class(time_scheme_t),       allocatable     :: time_scheme
        class(nonlinear_solver_t),  allocatable     :: nonlinear_solver
        class(linear_solver_t),     allocatable     :: linear_solver
        class(preconditioner_t),    allocatable     :: preconditioner

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
    !!  @param[in]  level   Initialization level specification. 'env', 'io', or 'finalize'
    !!
    !--------------------------------------------------------------------------------------------
    subroutine init(self,level)
        class(chidg_t),  intent(inout)   :: self
        character(*),    intent(in)      :: level


        ! Valid strings are:
        !   - 'env'         Basic environment initialization. Equations and supporting grid data
        !   - 'io'          Read namelist input.
        !   - 'finalize'    Call component initialization routines before run

        ! Call environment initialization routines by default on first init call
        if (.not. self%envInitialized ) then
            call log_init()


            ! Order matters here. Functions need to come first. Used by equations and bcs.
            call register_functions()
            call register_equations()
            call register_bcs()

            call initialize_grid()



            self%envInitialized = .true.
        end if





        select case (trim(level))

            !
            ! Do nothing else. Just to ensure the environment commands above are executed
            !
            case ('env')


            case ('mpi')
                call chidg_mpi_init()

            !
            ! Read Namelist input file
            !
            case ('io')
                call read_input()


            case ('communication')
                call establish_communication(self%data%mesh)

            !
            ! Initialize chimera
            !
            case ('chimera')
                call detect_chimera_faces(self%data%mesh)
                call detect_chimera_donors(self%data%mesh)
                call compute_chimera_interpolators(self%data%mesh)


            !
            ! Allocate components, based on input or default input data
            !
            case ('finalize')

                !
                ! Test chidg necessary components have been allocated
                !
                if (.not. allocated(self%time_scheme))      call chidg_signal(FATAL,"chidg%time_scheme component was not allocated")
                if (.not. allocated(self%nonlinear_solver)) call chidg_signal(FATAL,"chidg%nonlinear_solver component was not allocated")
                if (.not. allocated(self%linear_solver))    call chidg_signal(FATAL,"chidg%linearsolver component was not allocated")
                if (.not. allocated(self%preconditioner))   call chidg_signal(FATAL,"chidg%preconditioner component was not allocated")


                call self%time_scheme%init(self%data)
                call self%preconditioner%init(self%data)



            case default
                call chidg_signal(WARN,'chidg_t: Invalid initialization string')

        end select


    end subroutine init
    !**********************************************************************************************************








    !>  Set ChiDG environment components
    !!
    !!      -   Set time-scheme
    !!      -   Set nonlinear solver
    !!      -   Set linear solver
    !!      -   Set preconditioner
    !!      -   Set number of allocated domains (default=1)
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !!  @param[in]      selector    Character string for selecting the chidg component for initialization
    !!  @param[in]      selection   Character string for specializing the component being initialized
    !!  @param[inout]   options     Dictionary for initialization options
    !!
    !--------------------------------------------------------------------------------------------------------
    subroutine set(self,selector,selection,options)
        class(chidg_t),         intent(inout)   :: self
        character(*),           intent(in)      :: selector
        character(*),           intent(in)      :: selection
        type(dict_t), optional, intent(inout)   :: options 

        integer(ik) :: ierr



        select case (trim(selector))
            !
            ! Allocation for time scheme
            !
            case ('Time','time','time_scheme','Time_Scheme','timescheme','TimeScheme')
                if (allocated(self%time_scheme)) then
                    deallocate(self%time_scheme)
                    call create_time_scheme(selection,self%time_scheme,options)
                else
                    call create_time_scheme(selection,self%time_scheme,options)
                end if



            !
            ! Allocation for nonlinear solver
            !
            case ('nonlinearsolver','NonlinearSolver','nonlinear_solver','Nonlinear_Solver')
                if (allocated(self%nonlinear_solver)) then
                    deallocate(self%nonlinear_solver)
                    call create_nonlinear_solver(selection,self%nonlinear_solver,options)
                else
                    call create_nonlinear_solver(selection,self%nonlinear_solver,options)
                end if




            !
            ! Allocation for linear solver
            !
            case ('linearsolver','LinearSolver','linear_solver','Linear_Solver')
                if (allocated(self%linear_solver)) then
                    deallocate(self%linear_solver)
                    call create_linear_solver(selection,self%linear_solver,options)
                else
                    call create_linear_solver(selection,self%linear_solver,options)
                end if


            !
            ! Allocation for preconditioner
            !
            case ('preconditioner','Preconditioner')
                if (allocated(self%preconditioner)) deallocate(self%preconditioner)

                call create_preconditioner(selection,self%preconditioner)



            case default
                call chidg_signal_one(FATAL,"chidg%set: component string was not recognized. Check spelling and that the component was registered as an option in the chidg%set routine",selector)


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
    !!  @param[in]  partition   Optional partition for partial read of grid data during parallel execution.
    !!
    !!  TODO: Generalize spacedim
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine read_grid(self,gridfile,spacedim,partition)
        class(chidg_t),     intent(inout)           :: self
        character(*),       intent(in)              :: gridfile
        integer(ik),        intent(in)              :: spacedim
        type(partition_t),  intent(in), optional    :: partition

        character(len=5),   dimension(1)    :: extensions
        character(len=:),   allocatable     :: extension
        type(meshdata_t),   allocatable     :: meshdata(:)
        integer                             :: iext, extloc, idom, ndomains


        !
        ! Get filename extension
        !
        extensions = ['.h5']
        extension = get_file_extension(gridfile, extensions)


        !
        ! Call grid reader based on file extension
        !
        if ( extension == '.h5' ) then
            if (present(partition)) then
                call read_grid_partition_hdf(gridfile,partition,meshdata)
            else
                call read_grid_hdf(gridfile,meshdata)
            end if
        else
            call chidg_signal(FATAL,"chidg%read_grid: grid file extension not recognized")
        end if




        !
        ! Add domains to ChiDG%data
        !
        ndomains = size(meshdata)
        do idom = 1,ndomains

            call self%data%add_domain(                              &
                                      trim(meshdata(idom)%name),    &
                                      meshdata(idom)%points,        &
                                      meshdata(idom)%connectivity,  &
                                      spacedim,                     &
                                      meshdata(idom)%nterms_c,      &
                                      meshdata(idom)%eqnset         &
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
    !subroutine read_boundaryconditions(self, gridfile)
    subroutine read_boundaryconditions(self, gridfile, partition)
        class(chidg_t),     intent(inout)           :: self
        character(*),       intent(in)              :: gridfile
        type(partition_t),  intent(in), optional    :: partition

        character(len=5),   dimension(1)    :: extensions
        character(len=:),   allocatable     :: extension, dname
        type(bcdata_t),     allocatable     :: bcdata(:)
        integer                             :: idom, ndomains, iface, ierr

        !
        ! Get filename extension
        !
        extensions = ['.h5']
        extension = get_file_extension(gridfile, extensions)


        !
        ! Call boundary condition reader based on file extension
        !
        if ( extension == '.h5' ) then
            if (present(partition)) then
                call read_boundaryconditions_partition_hdf(gridfile,bcdata,partition)
            else
                call read_boundaryconditions_hdf(gridfile,bcdata)
            end if
        else
            call chidg_signal(FATAL,"chidg%read_boundaryconditions: grid file extension not recognized")
        end if



        !
        ! Add boundary conditions to ChiDG
        !
        ndomains = size(bcdata)
        do idom = 1,ndomains


            dname = bcdata(idom)%domain_
            do iface = 1,NFACES

                if ( allocated(bcdata(idom)%bcs(iface)%bc) ) then
                    if ( trim(bcdata(idom)%bcs(iface)%bc%name) /= "empty" ) then
                        call self%data%add_bc(dname, bcdata(idom)%bcs(iface)%bc, bcdata(idom)%bc_connectivity(iface))
                    end if
                end if

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
        class(chidg_t),     intent(inout)   :: self
        character(*),       intent(in)      :: solutionfile

        character(len=5),   dimension(1)    :: extensions
        character(len=:),   allocatable     :: extension
        type(meshdata_t),   allocatable     :: solutiondata(:)
        integer                             :: iext, extloc, idom, ndomains


        !
        ! Get filename extension
        !
        extensions = ['.h5']
        extension = get_file_extension(solutionfile, extensions)


        !
        ! Call grid reader based on file extension
        !
        if ( extension == '.h5' ) then
            call read_solution_hdf(solutionfile,self%data)
        else
            call chidg_signal(FATAL,"chidg%read_solution: grid file extension not recognized")
        end if


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
    subroutine initialize_solution_domains(self,nterms_s)
        class(chidg_t),     intent(inout)   :: self
        integer(ik),        intent(in)      :: nterms_s

        !
        ! TODO: put in checks for prerequisites
        !



        !
        ! Call domain solution storage initialization: mesh data structures that depend on solution expansion etc.
        !
        call self%data%initialize_solution_domains(nterms_s)


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
        class(chidg_t),     intent(inout)   :: self
        character(*),       intent(in)      :: solutionfile

        character(len=5),   dimension(1)    :: extensions
        character(len=:),   allocatable     :: extension
        type(meshdata_t),   allocatable     :: solutiondata(:)
        integer                             :: iext, extloc, idom, ndomains


        !
        ! Get filename extension
        !
        extensions = ['.h5']
        extension = get_file_extension(solutionfile, extensions)


        !
        ! Call grid reader based on file extension
        !
        if ( extension == '.h5' ) then
            call write_solution_hdf(solutionfile,self%data)
        else
            call chidg_signal(FATAL,"chidg%write_solution: grid file extension not recognized")
        end if


    end subroutine write_solution
    !************************************************************************************************************














    !>  Run ChiDG simulation
    !!
    !!      - This routine passes the domain data, nonlinear solver, linear solver, and preconditioner
    !!        components to the time scheme for iteration
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine run(self)
        class(chidg_t),     intent(inout)   :: self


        call self%time_scheme%iterate(self%data,self%nonlinear_solver,self%linear_solver,self%preconditioner)


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
    subroutine report(self)
        class(chidg_t), intent(inout)   :: self


        !call self%time_scheme%report()
        call self%nonlinear_solver%report()
        !call self%preconditioner%report()


    end subroutine report
    !************************************************************************************************************









    !> Any activities that need performed before the program completely terminates.
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
