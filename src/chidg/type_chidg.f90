module type_chidg
#include <messenger.h>
    use mod_equations,      only: initialize_equations
    use mod_grid,           only: initialize_grid
    use mod_io,             only: read_input

    use type_domain,        only: domain_t
    use atype_time_scheme,  only: time_scheme_t
    use atype_matrixsolver, only: matrixsolver_t
    use type_dict,          only: dict_t

    use mod_time_scheme,    only: create_time_scheme
    use mod_matrixsolver,   only: create_matrixsolver

    implicit none





    !>  The ChiDG Environment container
    !!
    !!      - Contains an array of domains, a time advancement scheme, and a matrix solver
    !!
    !!  @author Nathan A. Wukie
    !!
    !-----------------------------------------------------------------
    type, public    :: chidg_t

        type(domain_t),         allocatable     :: domains(:)
        class(time_scheme_t),   allocatable     :: timescheme
        class(matrixsolver_t),  allocatable     :: matrixsolver

        logical :: envInitialized = .false.
    contains
        procedure   :: init
        procedure   :: set
        procedure   :: run
    end type



contains

    !> ChiDG environment initialization routine
    !!      - Call initiailization procedures for equations, grid data, reading input
    !!
    !!  @author Nathan A. Wukie
    !!
    !!  @param[in]  level   Initialization level specification. 'env' or 'full'
    !---------------------------------------------------------------
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
                call self%timescheme%init(self%domains(1))




            case default
                call signal(WARN,'chidg_t: Invalid initialization string')

        end select


    end subroutine






    !>  Set ChiDG environment components
    !!
    !!      -   Set time-scheme
    !!      -   Set matrix solver
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
            case ('Time','time','time_scheme','Time_Scheme')
                if (allocated(self%timescheme)) then
                    deallocate(self%timescheme)
                    call create_time_scheme(selection,self%timescheme,options)
                else
                    call create_time_scheme(selection,self%timescheme,options)
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
            ! Allocation for number of domains
            !
            case ('ndomains','domains')
                if (allocated(self%domains)) then
                    call signal(WARN,"chidg: Domains already allocated")
                else
                    allocate(self%domains(1), stat=ierr)
                    if (ierr /= 0) call AllocationError
                end if




        end select


    end subroutine











    !>  Run ChiDG simultion
    !!
    !!      - This routine passes the domain and matrixsolver components to the
    !!        time scheme for iteration
    !!
    !!  @author Nathan A. Wukie
    !!
    !------------------------------------------------------------------------
    subroutine run(self)
        class(chidg_t),     intent(inout)   :: self


        call self%timescheme%solve(self%domains(1),self%matrixsolver)


    end subroutine












end module type_chidg
