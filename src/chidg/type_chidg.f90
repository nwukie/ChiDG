module type_chidg
#include <messenger.h>
    use mod_equations,          only: initialize_equations
    use mod_grid,               only: initialize_grid
    use mod_io,                 only: read_input

    use type_chidg_data,        only: chidg_data_t
    use type_timescheme,        only: timescheme_t
    use atype_matrixsolver,     only: matrixsolver_t
    use type_preconditioner,    only: preconditioner_t
    use type_dict,              only: dict_t

    use mod_timescheme,         only: create_timescheme
    use mod_matrixsolver,       only: create_matrixsolver
    use mod_preconditioner,     only: create_preconditioner

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
        procedure   :: set
        procedure   :: run
        procedure   :: report
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
            ! Allocate components, based on input or default input data
            !
            case ('finalize')

                ! Test chidg necessary components have been allocated
                if (.not. allocated(self%timescheme))     call signal(FATAL,"chidg%timescheme component was not allocated")
                if (.not. allocated(self%matrixsolver))   call signal(FATAL,"chidg%matrixsolver component was not allocated")
                if (.not. allocated(self%preconditioner)) call signal(FATAL,"chidg%preconditioner component was not allocated")


                call self%timescheme%init(self%data)
                call self%preconditioner%init(self%data)



            case default
                call signal(WARN,'chidg_t: Invalid initialization string')

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



            !
            ! Allocation for number of domains
            !
            !case ('ndomains','domains')
            !    if (allocated(self%domains)) then
            !        call signal(WARN,"chidg%set: Domains already allocated")
            !    else
            !        allocate(self%domains(1), stat=ierr)
            !        if (ierr /= 0) call AllocationError
            !    end if



            case default
                call signal(FATAL,"chidg%set: component string was not recognized. Check spelling and that the component was registered as an option in the chidg%set routine")


        end select


    end subroutine











    !>  Run ChiDG simulation
    !!
    !!      - This routine passes the domain and matrixsolver components to the
    !!        time scheme for iteration
    !!
    !!  @author Nathan A. Wukie
    !!
    !----------------------------------------------------------------------------
    subroutine run(self)
        class(chidg_t),     intent(inout)   :: self

        call self%timescheme%solve(self%data,self%matrixsolver,self%preconditioner)

    end subroutine run










    !> Report on ChiDG simulation
    !!
    !!  @author Nathan A. Wukie
    !!
    !!
    !!
    !!
    !!
    !----------------------------------------------------------------------------
    subroutine report(self)
        class(chidg_t), intent(inout)   :: self


        call self%timescheme%report()
        !call self%preconditioner%report()


    end subroutine report






end module type_chidg
