!> Chimera-based, discontinuous Galerkin equation solver
!!
!! This program is designed to solve partial differential equations,
!! and systems of partial differential equations, using the discontinuous
!! Galerkin method for spatial discretization using Chimera, overset grids to
!! represent the simulation domain.
!!
!! @author Nathan A. Wukie
!!
!!
!---------------------------------------------------------------------------------------------



program driver
    use mod_kinds,              only: rk, ik
    use mod_constants,          only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, ZERO
    use type_chidg,             only: chidg_t
    use mod_grid_operators,     only: initialize_variable
    use atype_function,         only: function_t
    use type_dict,              only: dict_t
    use mod_function,           only: create_function
    use mod_tecio,              only: write_tecio_variables
    use mod_io
    
    !
    ! Variable declarations
    !
    implicit none
    type(chidg_t)                       :: chidg
    class(function_t),  allocatable     :: constant, vortex, sod, roe
    type(dict_t)                        :: toptions, moptions



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
    call chidg%set('time_scheme',timescheme,toptions)
    call chidg%set('matrixsolver',matrixsolver,moptions)
    call chidg%set('preconditioner',preconditioner)




    !
    ! Add boundary conditions
    !
    call chidg%data%add_bc('D_01','euler_totalinlet',XI_MIN)
    call chidg%data%add_bc('D_01','euler_pressureoutlet',XI_MAX)
    call chidg%data%add_bc('D_01','euler_wall',ETA_MIN)
    call chidg%data%add_bc('D_01','euler_wall',ETA_MAX)
    call chidg%data%add_bc('D_01','euler_wall',ZETA_MIN)
    call chidg%data%add_bc('D_01','euler_wall',ZETA_MAX)




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
        call constant%set('val',1.13_rk)
        call initialize_variable(chidg%data,1,constant)

        ! rho_u
        call constant%set('val',190._rk)
        call initialize_variable(chidg%data,2,constant)

        ! rho_v
        call constant%set('val',ZERO)
        call initialize_variable(chidg%data,3,constant)

        ! rho_w
        call constant%set('val',ZERO)
        call initialize_variable(chidg%data,4,constant)

        ! rho_E
        call constant%set('val',248000._rk)
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





end program driver
