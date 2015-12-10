module type_chidg
#include <messenger.h>
    use mod_equations,          only: initialize_equations
    use mod_grid,               only: initialize_grid
    use mod_io,                 only: read_input, nterms_s, eqnset
    use mod_string_utilities,   only: get_file_extension

    use type_chidg_data,        only: chidg_data_t
    use type_timescheme,        only: timescheme_t
    use atype_matrixsolver,     only: matrixsolver_t
    use type_preconditioner,    only: preconditioner_t
    use type_meshdata,          only: meshdata_t
    use type_dict,              only: dict_t

    use mod_timescheme,         only: create_timescheme
    use mod_matrixsolver,       only: create_matrixsolver
    use mod_preconditioner,     only: create_preconditioner
    use mod_chimera,            only: detect_chimera_faces,  &
                                      detect_chimera_donors, &
                                      compute_chimera_interpolators

    use mod_hdfio,              only: read_grid_hdf, read_solution_hdf, write_solution_hdf

    implicit none





    !>  The ChiDG Environment container
    !!
    !!      - Contains an array of domains, a time advancement scheme, a matrix solver, and a preconditioner
    !!
    !!  @author Nathan A. Wukie
    !!
    !-----------------------------------------------------------------
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
        procedure   :: read_solution
        procedure   :: write_solution
    end type



contains

    !> ChiDG environment initialization routine
    !!      - Call initiailization procedures for equations, grid data, reading input
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  level   Initialization level specification. 'env', 'io', or 'finalize'
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
            call initialize_equations()
            call initialize_grid()
            call log_init()
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


            case ('chimera')
                call detect_chimera_faces(self%data%mesh)
                call detect_chimera_donors(self%data%mesh)
                call compute_chimera_interpolators(self%data%mesh)


            !
            ! Allocate components, based on input or default input data
            !
            case ('finalize')

                ! Test chidg necessary components have been allocated
                if (.not. allocated(self%timescheme))     call chidg_signal(FATAL,"chidg%timescheme component was not allocated")
                if (.not. allocated(self%matrixsolver))   call chidg_signal(FATAL,"chidg%matrixsolver component was not allocated")
                if (.not. allocated(self%preconditioner)) call chidg_signal(FATAL,"chidg%preconditioner component was not allocated")


                call self%timescheme%init(self%data)
                call self%preconditioner%init(self%data)



            case default
                call chidg_signal(WARN,'chidg_t: Invalid initialization string')

        end select


    end subroutine
    !------------------------------------------------------------------------------------------


















    !>  Set ChiDG environment components
    !!
    !!      -   Set time-scheme
    !!      -   Set matrix solver
    !!      -   Set preconditioner
    !!      -   Set number of allocated domains (default=1)
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]      selector    Character string for selecting the chidg component for initialization
    !!  @param[in]      selection   Character string for specializing the component being initialized
    !!  @param[inout]   options     Dictionary for initialization options
    !!
    !-----------------------------------------------------------------------
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
                call chidg_signal(FATAL,"chidg%set: component string was not recognized.      &
                                         Check spelling and that the component was registered &
                                         as an option in the chidg%set routine")


        end select


    end subroutine








    !>
    !!
    !!
    !!
    !!
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
                                      eqnset,                       &
                                      nterms_s                      &
                                      )
        end do




    end subroutine read_grid
    !##############################################################################################





    !>
    !!
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------
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
    !##############################################################################################






    !>
    !!
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------------------------
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
    !##############################################################################################









    !>  Run ChiDG simulation
    !!
    !!      - This routine passes the domain data, matrixsolver, and preconditioner
    !!        components to the time scheme for iteration
    !!
    !!  @author Nathan A. Wukie
    !!
    !---------------------------------------------------------------------------------------------
    subroutine run(self)
        class(chidg_t),     intent(inout)   :: self


        call self%timescheme%solve(self%data,self%matrixsolver,self%preconditioner)


    end subroutine run
    !##############################################################################################










    !> Report on ChiDG simulation
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine report(self)
        class(chidg_t), intent(inout)   :: self


        call self%timescheme%report()
        !call self%preconditioner%report()


    end subroutine report
    !################################################################################################









    !> Any activities that need performed before the program completely terminates.
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !-------------------------------------------------------------------------------------------------
    subroutine close(self)
        class(chidg_t), intent(inout)   :: self


        call log_finalize()


    end subroutine
    !##################################################################################################











end module type_chidg
