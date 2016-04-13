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
    use mod_kinds,              only: rk, ik
    use type_chidg,             only: chidg_t
    use mod_grid_operators,     only: initialize_variable
    use type_dict,              only: dict_t
    use type_function,          only: function_t
    use mod_function,           only: create_function
    use mod_tecio,              only: write_tecio_variables
    use mod_chidg_edit,         only: chidg_edit
    use mod_chidg_convert,      only: chidg_convert
    use mod_chidg_interpolate,  only: chidg_interpolate
    use mod_kirchoffs,          only: kirchoff
    use mod_io
    
    !
    ! Variable declarations
    !
    implicit none
    type(chidg_t)                       :: chidg
    type(dict_t)                        :: toptions, moptions
    class(function_t),  allocatable     :: constant

    integer(ik)                         :: narg
    character(len=1024)                 :: chidg_action, filename, file_a, file_b



    !
    ! Check for command-line arguments
    !
    narg = command_argument_count()


    !
    ! Execute ChiDG calculation
    !
    if ( narg == 0 ) then


        !
        ! Initialize ChiDG environment
        !
        call chidg%init('env')
        call chidg%init('io')


        !
        ! Read grid data from file
        !
        call chidg%read_grid(gridfile, spacedim)


        !
        ! Read boundary conditions
        !
        call chidg%read_boundaryconditions(gridfile)



        !
        ! Set time-scheme options
        !
        call toptions%set('dt',dt)
        call toptions%set('tol',ttol)
        call toptions%set('nsteps',nsteps)
        call toptions%set('nwrite',nwrite)
        call toptions%set('cfl0',cfl0)


        !
        ! Set matrix solver options
        !
        call moptions%set('tol',mtol)


        !
        ! Set ChiDG components
        !
        call chidg%set('timescheme',    timescheme,   toptions)
        call chidg%set('matrixsolver',  matrixsolver, moptions)
        call chidg%set('preconditioner',preconditioner)


        !
        ! Initialize solution data storage
        !
        call chidg%initialize_solution_domains(nterms_s)
        call chidg%init('chimera')



        !call chidg%data%init_sdata()
        call chidg%initialize_solution_solver()



        !
        ! Initialize solution
        !
        if (solutionfile_in == 'none') then
            call create_function(constant,'constant')


            ! rho
            call constant%set_option('val',0._rk)
            call initialize_variable(chidg%data,1,constant)

            ! rho_u
            call constant%set_option('val',0._rk)
            call initialize_variable(chidg%data,2,constant)

            ! rho_v
            call constant%set_option('val',0._rk)
            call initialize_variable(chidg%data,3,constant)

            ! rho_w
            call constant%set_option('val',0._rk)
            call initialize_variable(chidg%data,4,constant)

            ! rho_E
            call constant%set_option('val',0._rk)
            call initialize_variable(chidg%data,5,constant)


            ! rho
            call constant%set_option('val',0._rk)
            call initialize_variable(chidg%data,6,constant)

            ! rho_u
            call constant%set_option('val',0._rk)
            call initialize_variable(chidg%data,7,constant)

            ! rho_v
            call constant%set_option('val',0._rk)
            call initialize_variable(chidg%data,8,constant)

            ! rho_w
            call constant%set_option('val',0._rk)
            call initialize_variable(chidg%data,9,constant)

            ! rho_E
            call constant%set_option('val',0._rk)
            call initialize_variable(chidg%data,10,constant)




        else

            !
            ! TODO: put in check that solutionfile actually contains solution
            !
            call chidg%read_solution(solutionfile_in)

        end if

        

        !
        ! Wrap-up initialization activities
        !
        call chidg%init('finalize')



        !
        ! Write initial solution
        !
        if (initial_write) call chidg%write_solution(solutionfile_out)





        !
        ! Run ChiDG simulation
        !
        call chidg%run()


        !
        ! Write final solution
        !
        if (final_write) call chidg%write_solution(solutionfile_out)




        !
        ! Reporting
        !
        call chidg%report()





        !
        ! Close ChiDG
        !
        call chidg%close()






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

        else if ( trim(chidg_action) == 'kirchoff' ) then
            call kirchoff(filename)


        else
            call chidg_signal(FATAL,"chidg: unrecognized action '"//trim(chidg_action)//"'. Valid options are: 'edit', 'convert'")

        end if



    !
    ! ChiDG tool execution. 3 arguments.
    !
    else if ( narg == 3 ) then


        call get_command_argument(1,chidg_action)
        call get_command_argument(2,file_a)
        call get_command_argument(3,file_b)
        

        if ( trim(chidg_action) == 'interpolate' ) then
            call chidg_interpolate(trim(file_a), trim(file_b))

        else
            call chidg_signal(FATAL,"chidg: unrecognized action '"//trim(chidg_action)//"'. Valid options are: 'edit', 'convert'")

        end if







    else
        call chidg_signal(FATAL,"chidg: invalid number of arguments. Expecting (0) arguments: 'chidg'. or (2) arguments: 'chidg action file'.")
    end if








end program driver
