module type_chidg
#include <messenger.h>
    use mod_constants,          only: NFACES
    use mod_equations,          only: register_equations
    use mod_bc,                 only: register_bcs
    use mod_function,           only: register_functions
    use mod_grid,               only: initialize_grid
    use mod_io,                 only: read_input, nterms_s
    use mod_string_utilities,   only: get_file_extension

    use type_chidg_data,        only: chidg_data_t
    use type_timescheme,        only: timescheme_t
    use type_matrixsolver,      only: matrixsolver_t
    use type_preconditioner,    only: preconditioner_t
    use type_meshdata,          only: meshdata_t
    use type_bcdata,            only: bcdata_t
    use type_dict,              only: dict_t

    use mod_timescheme,         only: create_timescheme
    use mod_matrixsolver,       only: create_matrixsolver
    use mod_preconditioner,     only: create_preconditioner
    use mod_chimera,            only: detect_chimera_faces,  &
                                      detect_chimera_donors, &
                                      compute_chimera_interpolators

    use mod_hdfio,              only: read_grid_hdf, read_boundaryconditions_hdf, read_solution_hdf, write_solution_hdf

    implicit none





    !>  The ChiDG Environment container
    !!
    !!      - Contains an array of domains, a time advancement scheme, a matrix solver, and a preconditioner
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !--------------------------------------------------------------------------------------------------------
    type, public    :: chidg_t

        type(chidg_data_t)                          :: data
        class(timescheme_t),        allocatable     :: timescheme
        class(matrixsolver_t),      allocatable     :: matrixsolver
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

            !
            ! Order matters here. Functions need to come first. Used by equations and bcs.
            !
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


            !
            ! Read Namelist input file
            !
            case ('io')
                call read_input()


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
                if (.not. allocated(self%timescheme))     call chidg_signal(FATAL,"chidg%timescheme component was not allocated")
                if (.not. allocated(self%matrixsolver))   call chidg_signal(FATAL,"chidg%matrixsolver component was not allocated")
                if (.not. allocated(self%preconditioner)) call chidg_signal(FATAL,"chidg%preconditioner component was not allocated")


                call self%timescheme%init(self%data)
                call self%preconditioner%init(self%data)



            case default
                call chidg_signal(WARN,'chidg_t: Invalid initialization string')

        end select


    end subroutine init
    !**********************************************************************************************************








    !>  Set ChiDG environment components
    !!
    !!      -   Set time-scheme
    !!      -   Set matrix solver
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
                if (allocated(self%timescheme)) then
                    deallocate(self%timescheme)
                    call create_timescheme(selection,self%timescheme,options)
                else
                    call create_timescheme(selection,self%timescheme,options)
                end if


            !
            ! Allocation for matrix solver
            !
            case ('matrixsolver','MatrixSolver','matrix_solver','Matrix_Solver')
                if (allocated(self%matrixsolver)) then
                    deallocate(self%matrixsolver)
                    call create_matrixsolver(selection,self%matrixsolver,options)
                else
                    call create_matrixsolver(selection,self%matrixsolver,options)
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
    !!
    !----------------------------------------------------------------------------------------------
    subroutine read_grid(self,gridfile)
        class(chidg_t),     intent(inout)   :: self
        character(*),       intent(in)      :: gridfile

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
            call read_grid_hdf(gridfile,meshdata)
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
                                      meshdata(idom)%nterms_c,      &
                                      meshdata(idom)%eqnset,        &
                                      nterms_s                      &
                                      )
        end do


    end subroutine read_grid
    !*************************************************************************************************************











    !>  Read boundary conditions from grid file.
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/5/2016
    !!
    !!  @param[in]  gridfile    String specifying a gridfile, including extension.
    !!
    !-------------------------------------------------------------------------------------------------------------
    subroutine read_boundaryconditions(self, gridfile)
        class(chidg_t), intent(inout)   :: self
        character(*),   intent(in)      :: gridfile

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
            call read_boundaryconditions_hdf(gridfile,bcdata)
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

                    call self%data%add_bc(dname, bcdata(idom)%bcs(iface)%bc, bcdata(idom)%bcface(iface))

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
    !!      - This routine passes the domain data, matrixsolver, and preconditioner
    !!        components to the time scheme for iteration
    !!
    !!  @author Nathan A. Wukie
    !!  @date   2/1/2016
    !!
    !------------------------------------------------------------------------------------------------------------
    subroutine run(self)
        class(chidg_t),     intent(inout)   :: self


        call self%timescheme%solve(self%data,self%matrixsolver,self%preconditioner)


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


        call self%timescheme%report()
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
    subroutine close(self)
        class(chidg_t), intent(inout)   :: self


        call log_finalize()


    end subroutine close
    !************************************************************************************************************











end module type_chidg
