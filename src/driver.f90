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
    use mod_constants,          only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, ZERO
    use type_chidg,             only: chidg_t
    use mod_grid_operators,     only: initialize_variable
    use type_dict,              only: dict_t
    use type_function,          only: function_t
    use mod_function,           only: create_function
    use type_bc,                only: bc_t
    use mod_bc,                 only: create_bc
    use mod_tecio,              only: write_tecio_variables
    use mod_io
    use mod_chidg_edit,         only: chidg_edit
    
    !
    ! Variable declarations
    !
    implicit none
    type(chidg_t)                       :: chidg
    type(dict_t)                        :: toptions, moptions

    class(bc_t),        allocatable     :: bc_wall, bc_inlet, bc_outlet
    class(function_t),  allocatable     :: constant, vortex, sod, roe

    integer(ik)                         :: narg
    !character(len=:),   allocatable     :: chidg_action, filename
    character(len=1024)                 :: chidg_action, filename




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
        call chidg%read_grid(gridfile)




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
        ! Set up boundary conditions
        !
        call create_bc('euler_wall',            bc_wall  )
        call create_bc('euler_totalinlet',      bc_inlet )
        call create_bc('euler_pressureoutlet',  bc_outlet)


        call bc_outlet%set_fcn(       'Static Pressure','constant')
        call bc_outlet%set_fcn_option('Static Pressure','val',107000._rk)



        !
        ! Add boundary conditions
        !
        !call chidg%data%add_bc('D_01','euler_totalinlet',XI_MIN)
        !call chidg%data%add_bc('D_01','euler_pressureoutlet',XI_MAX)
        call chidg%data%add_bc('D_01',bc_wall,ETA_MIN)
        call chidg%data%add_bc('D_01',bc_wall,ETA_MAX)
        call chidg%data%add_bc('D_01',bc_wall,ZETA_MIN)
        call chidg%data%add_bc('D_01',bc_wall,ZETA_MAX)


        !call chidg%data%add_bc('D_02','euler_totalinlet',XI_MIN)
        !call chidg%data%add_bc('D_02','euler_pressureoutlet',XI_MAX)
        call chidg%data%add_bc('D_02',bc_outlet,ETA_MIN)
        call chidg%data%add_bc('D_02',bc_wall,ETA_MAX)
        call chidg%data%add_bc('D_02',bc_wall,ZETA_MIN)
        call chidg%data%add_bc('D_02',bc_wall,ZETA_MAX)


       !call chidg%data%add_bc('D_03','euler_totalinlet',XI_MIN)
       !call chidg%data%add_bc('D_03','euler_pressureoutlet',XI_MAX)
       call chidg%data%add_bc('D_03',bc_wall,ETA_MIN)
       call chidg%data%add_bc('D_03',bc_wall,ETA_MAX)
       call chidg%data%add_bc('D_03',bc_wall,ZETA_MIN)
       call chidg%data%add_bc('D_03',bc_wall,ZETA_MAX)


       !call chidg%data%add_bc('D_04','euler_totalinlet',XI_MIN)
       !call chidg%data%add_bc('D_04','euler_pressureoutlet',XI_MAX)
       call chidg%data%add_bc('D_04',bc_inlet,ETA_MIN)
       call chidg%data%add_bc('D_04',bc_wall,ETA_MAX)
       call chidg%data%add_bc('D_04',bc_wall,ZETA_MIN)
       call chidg%data%add_bc('D_04',bc_wall,ZETA_MAX)










        !
        ! Initialize solution data storage
        !
        call chidg%init('chimera')
        call chidg%data%init_sdata()




        !
        ! Initialize solution
        !
        if (solutionfile_in == 'none') then
            call create_function(constant,'constant')


            ! rho
            !call constant%set_option('val',1.13_rk)
            call constant%set_option('val',1.25_rk)
            call initialize_variable(chidg%data,1,constant)

            ! rho_u
            !call constant%set_option('val',190._rk)
            call constant%set_option('val',80._rk)
            call initialize_variable(chidg%data,2,constant)

            ! rho_v
            call constant%set_option('val',ZERO)
            call initialize_variable(chidg%data,3,constant)

            ! rho_w
            call constant%set_option('val',ZERO)
            call initialize_variable(chidg%data,4,constant)

            ! rho_E
            !call constant%set_option('val',248000._rk)
            call constant%set_option('val',270000._rk)
            call initialize_variable(chidg%data,5,constant)

        else

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
    ! ChiDG tool execution
    !
    else if ( narg == 2 ) then




        call get_command_argument(1,chidg_action)
        call get_command_argument(2,filename)
        chidg_action = trim(chidg_action)
        filename = trim(filename)
        

        if ( trim(chidg_action) == 'edit' ) then
            
            call chidg_edit(trim(filename))
        else
            call chidg_signal_one(FATAL,"chidg: unrecognized action.",trim(chidg_action))
        end if


    else
        call chidg_signal(FATAL,"chidg: invalid number of arguments. Expecting (0) arguments: 'chidg'. or (2) arguments: 'chidg action file'.")
    end if



































end program driver
