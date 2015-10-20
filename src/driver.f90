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
    use mod_constants,          only: XI_MIN, XI_MAX, ETA_MIN, ETA_MAX, ZETA_MIN, ZETA_MAX, ONE, TWO, THREE, FOUR, FIVE, ZERO
    use type_chidg,             only: chidg_t
    use type_meshdata,          only: meshdata_t
    use mod_hdfio,              only: read_grid_hdf
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
    type(meshdata_t),       allocatable :: meshdata(:)
    class(function_t),      allocatable :: constant, vortex, sod, roe
    type(dict_t)                        :: toptions
    integer(ik)                         :: nterms_c, idom, ndomains



    !
    ! Initialize ChiDG environment
    !
    call chidg%init('env')
    call chidg%init('io')

    print*, 'driver - 1'

    !
    ! Read grid data from file
    !
    call read_grid_hdf(gridfile,meshdata)

    print*, 'driver - 2'

    !
    ! Add domains to ChiDG
    !
    ndomains = size(meshdata)
    do idom = 1,ndomains
        call chidg%data%add_domain(trim(meshdata(idom)%name),meshdata(idom)%points,meshdata(idom)%nterms_c,eqnset,nterms_s)
    end do
    print*, 'driver - 3'

    !
    ! Initialize solution data storage
    !
    call chidg%data%init_sdata()


    print*, 'driver - 4'
    !
    ! Set time-scheme options
    !
    call toptions%set('dt',dt)
    call toptions%set('tol',ttol)
    call toptions%set('nsteps',nsteps)
    call toptions%set('nwrite',nwrite)
    call toptions%set('cfl0',cfl0)

    print*, 'driver - 5'
    !
    ! Set ChiDG components
    !
    call chidg%set('time_scheme',timescheme,toptions)
    call chidg%set('matrixsolver',matrixsolver)
    call chidg%set('preconditioner',preconditioner)

    print*, 'driver - 6'


    do idom = 1,chidg%data%ndomains
    associate ( data => chidg%data )
        !
        ! Initialize domain
        !
        call data%add_bc('D_01','euler_totalinlet',XI_MIN)
        call data%add_bc('D_01','euler_pressureoutlet',XI_MAX)
        call data%add_bc('D_01','euler_wall',ETA_MIN)
        call data%add_bc('D_01','euler_wall',ETA_MAX)
        call data%add_bc('D_01','euler_wall',ZETA_MIN)
        call data%add_bc('D_01','euler_wall',ZETA_MAX)

    print*, 'driver - 7'

        !
        ! Initialize solution
        !
        call create_function(constant,'constant')
        call create_function(sod,'sod')
        call create_function(roe,'roe_check') 

    
    print*, 'driver - 8'
        ! rho
        call constant%set('val',1.13262_rk)
        call initialize_variable(chidg%data,1,constant)

        ! rho_u
        call constant%set('val',190.339029_rk)
        call initialize_variable(chidg%data,2,constant)

        ! rho_v
        call constant%set('val',ZERO)
        call initialize_variable(chidg%data,3,constant)

        ! rho_w
        call constant%set('val',ZERO)
        call initialize_variable(chidg%data,4,constant)

        ! rho_E
        call constant%set('val',248493.425_rk)
        call initialize_variable(chidg%data,5,constant)
    print*, 'driver - 9'

    end associate
    end do 

    print*, 'driver - 10'
    !
    ! Write initial solution
    !
    if (initial_write) then
        call write_tecio_variables(chidg%data,'0.plt',1)
    end if



    print*, 'driver - 11'
    

    !
    ! Wrap-up initialization activities
    !
    call chidg%init('finalize')

    print*, 'driver - 12'




    !
    ! Run ChiDG simulation
    !
    call chidg%run()


    print*, 'driver - 13'



    !
    ! Reporting
    !
    call chidg%report()


    print*, 'driver - 14'

    !
    ! Write final solution
    !
    if (final_write) then
        call write_tecio_variables(chidg%data,'9999999.plt',1)
    end if


end program driver
